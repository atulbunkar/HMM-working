// HMM.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "vector"
#include "iostream"
#include <stdlib.h> 

#include <conio.h>
#include <sstream>
#include <list>

#pragma warning(disable:4996)

using namespace std;
#define ld long double


#define T 60 			//observations per sequence / frames to consider

#define N 5				//5 states
#define M 32			//32 possible outcomes
#define p 12			//12 CC
#define F 320			//frame size
#define K 32			// k of k means
#define slide 80		//how much to slide the window

#define unisize 10*20*T	  //Total frames/universe

double norm_data[45000];	//store data sample for one utterance.

double Uni[10][20+1][T][p+1] = {0} ;	//store final Ci for all digits. Feed to LBG
double Uni2d[unisize][p] = {0} ;	    //flattening Uni to 2D.
double test_Ci[T][p] = {0};				//storing test ci;

double codebook[K][p] = {0} ;			//store codebook

//for 10 files , 20x85 Ci 
int obs[unisize] = {0};
int O[T+1] = {0} ;	  //store the observation seq of all digit - 1 observation for 1 Ci .
int obs2d[20*10][T] ;

ld Pi[N+1];				//store pi
ld A[N+1][N+1];			//store A, state transition mtrx
ld B[N+1][M+1];			//store B, output probability mtrx

ld A_bar[N+1][N+1];			//store fixed/improved model
ld B_bar[N+1][M+1];		
ld Pi_bar[N+1];

ld Pi_avg[N+1];				//To store avg of 20 models
ld A_avg[N+1][N+1];			
ld B_avg[N+1][M+1];

//ld P ;					//store prob of a word belonging to a model
ld p_star = 0;					  //store p* 		
ld prev_p_star =0;
int q_star_t[T+1] = {0};      //store q_t* 

ld alpha[T+1][N+1] ; //forward proc
ld beta[T+1][N+1] ;  //backward proc
ld gamma[T+1][N+1] ; //gamma proc
ld delta[T+1][N+1] ; //viterbi
int psi[T+1][N+1] ;	 //viterbi
ld zeta[T+1][T+1][N+1] ; //baum welsh
int check = 0;

#include "utility.h"
#include "get_ceps.h"
#include "lbg.h"

//step 1
void get_Uni(){		//Get Ci of all digit samples
	char filenum[3] ;	
	int k =0 , i = 0 ,j =0 ; 

	for(int d = 0 ; d<=9 ; d++){
		for(int u = 1 ; u <= 20; u++){
			char path[30] = "./digits/";    
			sprintf(filenum ,"%d", d) ;
			strcat(path, filenum) ;
			strcat(path, "_") ;
			sprintf(filenum ,"%d", u) ;
			strcat(path, filenum) ;
			strcat(path, ".txt") ;
			//cout << path << endl; 
			genrate_Uni(d , u , path);
		}
	}
}


void dump_Uni(){		//dump universe to csv. Only for debugging check
	FILE* fp_dump;
	fp_dump = fopen("Universe.csv", "w+") ;
	for(int d = 0 ; d<=9 ; d++){
		for(int u = 1 ; u <= 20; u++){
			for(int frame = 0; frame < T ; frame++){
				for(int x=1; x <= p ; x++){
					fprintf(fp_dump,"%lf,",Uni[d][u][frame][x]) ;
				}
				fprintf(fp_dump , "\n");
			}
		}
	}
	fclose(fp_dump);
	return;
}

void Uni_2d(){						// flatten it b4 calling LBG.
	int row = 0 ;
	for(int d = 0 ; d<=9 ; d++){
		for(int u = 1 ; u <= 20; u++){
			for(int frame = 0; frame < T ; frame++){
				for(int x=1; x <= p ; x++){
					Uni2d[row][x-1] = Uni[d][u][frame][x];
				}
				row++;
			}
		}
	}
	return ;
}

void gen_cb(){	
	get_codebook(); //get cb from lbg (not used in assignment)
}

void dump_cb(){	//dump created codebook (not used in assignment)
	FILE* fp_dump ;
	fp_dump = fopen("Codebook.csv", "w+");
	for(int k = 0 ; k <32 ; k++){
		for(int i = 0; i <p ; i++){
			fprintf(fp_dump,"%lf,",codebook[k][i]) ;
		}
		fprintf(fp_dump , "\n");
	}
	fclose(fp_dump);	
	return ;
}

void read_cb(){	//read created codebook csv (not used)
	FILE* fp ;
	fp = fopen("Codebook.csv", "r");
	char line[200];
	double temp = 0;
	int row = 0, col = 0;

	while(fgets(line , sizeof(line) , fp)){  // reading CSV row by row
		char* tok ;
		tok = strtok(line , ",");			// tokenzing by delimiter 
		
		while(tok != NULL){	
			temp = atof(tok);	
			codebook[row][col] = temp;			//Reading Universe to 2d file.
			tok = strtok(NULL, ",");
			col++;
		}
		col = 0;
		row++;
	}
	fclose(fp);
	return ;
}

void read_cbtxt(){	// Codebook to use from TA (Used this one)
	string input; double number;

	char fullname[100] = "codebook.txt" ; 

    ifstream file(fullname);
	if(!file)cout<< "\n cant open ";

	int idx = 0 ; int row = 0;  

	while(getline(file, input)){
		istringstream iss( input );
		idx = 0;

		while( iss >> number ){
			 codebook[row][idx++] =number;
			 //cout<< codebook[row][idx-1] << " " ;
		}
		row++;
	}
	//cout << row;
	file.close();
	return ;

}

void get_obseravtion(){		//finds min dist index from codebook for all Universe - Creates Observation sequence
	double tok_wt[p] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

	for(int i = 0 ; i < unisize ; i++){			//For each of Uni.

		double minDist = DBL_MAX ; int minDistIndx = -1 ;
		for(int j = 0; j< K; j++){					//iterate each of CB.
			double temp = 0;

			for(int ci = 0; ci < p ; ci++){			//cur frame

				temp += ( tok_wt[ci] * pow((Uni2d[i][ci] - codebook[j][ci] ),2) ) ;
			}
			if(temp < minDist){
				minDist=temp;
				minDistIndx=j+1;
			}
		}
		obs[i] = minDistIndx;
		//cout<< i+1 << " " << obs[i] << endl;
	}
	return ;
}

void get_O(int digit ,int utter){		//gets Obs seq for one utterance from obs after filling obs from get_observation()
	
	int ind = 1;						//digit - > 0-9 , utter -> 1-20.
	int s = T*(20*digit + utter-1) ;			

	for(int i = s ; i < s+T ; i++){
		O[ind++] = obs[i];								//store O[1] to O[T] 
	}
	return;
}

void dump_obs(){	//dump to file
	FILE* fp_dump ;
	fp_dump = fopen("Obs Seq.csv", "w+");
	int idx = 0 ;

	for(int d = 0 ; d<=9 ; d++){				//prints 200 lines , each of T obs seq.
		for(int u = 1 ; u <= 20; u++){
			for(int frame = 0; frame < T ; frame++){
					fprintf(fp_dump,"%d,", obs[idx++]) ;
			}
			fprintf(fp_dump , "\n");
		}
	}
	fclose(fp_dump);	
	return;
}

void read_obs(){ //read obs seq to obs2d
	FILE* fp ;
	fp = fopen("Obs Seq.csv", "r");
	char line[200];
	int temp = 0;
	int row = 0,col = 0;

	while(fgets(line , sizeof(line) , fp)){  // reading CSV row by row
		char* tok ;
		tok = strtok(line , ",");			// tokenzing by delimiter 
		
		while(tok != NULL){	
			temp = atoi(tok);	
			obs2d[row][col] = temp;			//Reading Universe to 2d file.
			tok = strtok(NULL, ",");
			col++;
		}
		col = 0; row++;
	}
	return ;
	fclose(fp);
}

void get_O2(int d, int u){		//fill O from obs2d got from reading obs seq.csv
	int idx = 1;
	int row = 20*d + u-1;

	for(int i = 0 ; i < T ; i++){
		O[idx++] = obs2d[row][i];						
	}
	return;
}



//---------------------------------------------------------------------------
//MODEL BUILD
void feed_forward_model() //this is biased model 
{
	for(int i=1;i<=N;i++) //assign the given values
		if(i==1) //make the first state as starting state
			Pi[i]=1.0;
		else
			Pi[i]=0;
    for(int i=1;i<=N;i++) 
        for(int j=1;j<=N;j++)
			if(i==j&&i!=N)
				A[i][j]=0.8; 
			else if(i==j&&i==N)
				A[i][j]=1;
			else if(j==i+1)
				A[i][j]=0.2; 
			else
				A[i][j]=0; 
    for(int i=1;i<=N;i++) 
        for(int j=1;j<=M;j++)
            B[i][j]=1.0/M;

	return ;
}

//calculating alpha
ld forward_proc() //alpha computation
{
    for(int i=1;i<=N;i++) //initialization
        alpha[1][i]=Pi[i]*B[i][O[1]];
    for(int t=1;t<=T-1;t++) 
    {
        for(int j=1;j<=N;j++)
        {
            ld sum=0;
            for(int i=1;i<=N;i++)
            {
                sum+=alpha[t][i]*A[i][j];
            }
            alpha[t+1][j]=sum*B[j][O[t+1]];
        }
    }
    ld P=0;
    for(int i=1;i<=N;i++) //estimate what is the probability that the observation sequence is from the current model
	{
        P+=alpha[T][i];
	}
    return P;
}

//calculating beta
void backward_proc() //beta computation
{
    for(int i=1;i<=N;i++) //initialization
        beta[T][i]=1;

    for(int t=T-1;t>=1;t--) //induction
    {
        for(int i=1;i<=N;i++)
        {
            ld sum=0;
            for(int j=1;j<=N;j++)
            {
                sum+=A[i][j]*B[j][O[t+1]]*beta[t+1][j];
            }
            beta[t][i]=sum;
        }
    }
    return;
}

//calculating gamma
void gamma_proc() //gamma computation
{
	int *q,argmax=1; 
	q=new int[T+1];
	ld devider=0; 
	for(int t=1;t<=T;t++)
	{
		for(int i=1;i<=N;i++) //compute it once for t
		{
			devider+=alpha[t][i]*beta[t][i];
		}
		argmax=1;
		for(int i=1;i<=N;i++)
		{
			gamma[t][i]=alpha[t][i]*beta[t][i]/devider;
			if(gamma[t][argmax]<gamma[t][i])
				argmax=i;
		}
		q[t]=argmax;
		devider=0;
	}

	ld P=1;
	for(int t=1;t<=T;t++)
	{
		P*=gamma[t][q[t]];
		//cout << gamma[t][q[t]] << endl;
	}
	//cout << P << " " ;
	return;
}

//Soln-2 -calulating delta and psi
ld viterbi(){
	memset(q_star_t, 0, sizeof q_star_t);
	int argmax=0;			 // to store the argument which gives maximum probability
	for(int i=1 ;i<=N;i++)  //step 1 - initialization
	{
		delta[1][i]=Pi[i]*B[i][O[1]];
		//cout  << delta[1][i] << endl;
		psi[1][i]=0; //indicates, no state is yet assigned
	}
		
	//step 2 - Recursion
	for(int t=2;t<=T;t++) //recursion over time sequence
	{
		for(int j=1;j<=N;j++) // to store maximum probabilities
		{
			argmax=1; //assume the first state gives maximum probability
			for(int i=1;i<=N;i++)
			{
				if(delta[t-1][i]*A[i][j] > delta[t-1][argmax]*A[argmax][j]) //checking for maximum score
					argmax=i;							//argument which gives maximum score
			}
			delta[t][j] = delta[t-1][argmax]*A[argmax][j]*B[j][O[t]];
			psi[t][j]= argmax;					  //state which gives maximum probability
			//cout  << delta[t][j] << "  ";
		}
	}

	//step 3 - Termination
	argmax=1 ;
	for(int i=1; i <= N; i++){                 //to find the argmax for last time frame
		if(delta[T][i] > delta[T][argmax])
			argmax=i;
	}
	p_star= delta[T][argmax];				//holds best possible state's value

	q_star_t[T] = argmax;					//holds best possible state at last
	//cout << ob << " ---> " << q_star_t[0][os][T] << endl;
	//step 4 - Back Tracking the path
	for(int t=T-1; t>=1 ;t--){
		q_star_t[t]= psi[t+1][q_star_t[t+1]];
		//cout << t+1 << " --> " << q_star_t[0][os][t+1] << endl;
	}
	return p_star;
}

//soln 3 - BaumWelsh - calculating Zeta
void baum_welch_proc() //zeta computaion
{
	ld devider=0; //used as common devider just like in baye's theorem
	for(int t=1;t<=T-1;t++) //repeat this for all T-1 state transitions
	{
		devider=0;
		for(int i=1;i<=N;i++)
		{
			for(int j=1;j<=N;j++)
				devider+=alpha[t][i]*A[i][j]*B[j][O[t+1]]*beta[t+1][j];  //Need to confirm if correct. 
		}
		for(int i=1;i<=N;i++)
		{
			for(int j=1;j<=N;j++)
				zeta[t][i][j]=(alpha[t][i]*A[i][j]*B[j][O[t+1]]*beta[t+1][j]) / devider;
		}
	}
}

void re_estimation() //re-estimate transition probabilities and observation symbol probabilities
{
	ld numerator = 0, denominator =0 ; //for re-estimation of transition probabilities
	for(int i=1;i<=N;i++) //re-estimation of Pi
		Pi_bar[i]=gamma[1][i];

	for(int i=1;i<=N;i++) //re-estimation of A
	{
		for(int j=1;j<=N;j++)
		{
			numerator=0;
			denominator=0;
			for(int t=1;t<=T-1;t++)	  
			{
				numerator+=zeta[t][i][j];
				denominator+=gamma[t][i];
			}
			A_bar[i][j]=numerator/denominator;
		}
	}
	ld maxValue=-1, sum=0;
	int  maxIndex=-1;
	for(int j=1;j<=N;j++) //re-estimation of B
	{
		maxValue=-1, maxIndex=-1, sum=0;
		for(int k=1;k<=M;k++)
		{
			numerator=0;
			denominator=0;
			for(int t=1;t<=T;t++)
			{
				if(O[t]==k)
					numerator+=gamma[t][j];

				denominator+=gamma[t][j];
			}
			if(denominator==0){
				//cout<<"ERROR : DENOM = 0 !!"<<endl; 
				B_bar[j][k] = 0;
			}
		
			else{
				B_bar[j][k]=numerator/denominator;
			}
			if(B_bar[j][k] > maxValue){
					maxValue = B_bar[j][k] ;
					maxIndex = k;
			}
			if(B_bar[j][k] < (ld)(pow((ld)10,(ld)-30))){		//< (ld)(pow((ld)10,(ld)-30)) ->replace with-> ==0
				sum += (ld)(pow((ld)10,(ld)-30)) - B_bar[j][k] ;
				B_bar[j][k] = (ld)(pow((ld)10,(ld)-30)) ;
			}
		}
		B_bar[j][maxIndex] -= sum  ;
	}
}

void replace_model() //replace old model with new model
{
	//cout << "replacing the old model with new model" << endl;
	for(int i=1;i <= N;i++) 
        Pi[i]=Pi_bar[i];
	
    for(int i=1;i <= N;i++) 
        for(int j=1;j <= N;j++)
            A[i][j]=A_bar[i][j];

    for(int i=1;i <= N;i++) 
        for(int j=1;j <= M;j++)
            B[i][j]=B_bar[i][j];

	return;
}

//adding all converged model to avg
void add2avg(){

	for(int i=1;i<=N;i++)  //display initial state state distributions
	{
        Pi_avg[i] += Pi[i] ;
	}
    for(int i=1;i<=N;i++) //display state transition probabilities
	{
        for(int j=1;j<=N;j++)
		{
            A_avg[i][j] += A[i][j] ;	
		}
	
	}
    for(int i=1;i<=N;i++) //display observation symbol probability distribution
	{
        for(int j=1;j<=M;j++)
		{
            B_avg[i][j] += B[i][j] ;
		}
	}
	return;
}

//taking average of all models
void do_avg(int count){
	ld sum = 0; ld diff;
	for(int i=1;i<=N;i++)  
	{
        Pi_avg[i] = int( Pi_avg[i] / count ) ;
	}
    for(int i=1;i<=N;i++) 
	{
		sum = 0;
		int  maxIndex=-1; ld maxValue=-1;

        for(int j=1;j<=N;j++)
		{
            A_avg[i][j] = A_avg[i][j]/count*1.0 ;
			sum += A_avg[i][j];	//row-wise sum.

			if(A_avg[i][j] > maxValue){
				maxValue = A_avg[i][j] ;
				maxIndex = j;
			}
		}
		if(sum != 1){	//if sum != 1 , we add the diff to largest value.
			diff = 1 - sum;
			A_avg[i][maxIndex] = A_avg[i][maxIndex] + diff;
		}
	
	}


    for(int i=1;i<=N;i++) 
	{
		sum = 0;
		int  maxIndex=-1; ld maxValue=-1;

        for(int j=1;j<=M;j++)
		{
            B_avg[i][j] = B_avg[i][j]/count*1.0 ;
			sum += B_avg[i][j];	

			if(B_avg[i][j] > maxValue){
				maxValue = B_avg[i][j] ;
				maxIndex = j;
			}
		}
		if(sum != 1){	//if sum != 1 , we add the diff to largest value.
			diff = 1 - sum ;
			B_avg[i][maxIndex] = B_avg[i][maxIndex] + diff;
		}
		//cout << sum << endl;
	}
	return;
}

//Main process. Converging the models for given iterations
void RunModel(int iters){
	int itr = 0;         // iteration no.
	p_star=0;
	prev_p_star=-1;

	while(itr<=iters){
		itr++;
		forward_proc();
		backward_proc();

		gamma_proc();
		viterbi() ;

		baum_welch_proc();

		re_estimation();		// get re estimated model for current iteration
		replace_model();		//store re estimated model as original model to work on them again.
		
		if(abs(p_star - prev_p_star ) < (ld)(pow((ld)10,(ld)-50)) )break;  //check
		prev_p_star = p_star ;
	}
	cout << "final p_star: " <<  p_star << endl;

	cout<< "Final State Sequence : " << endl;
	for(int x = 1 ; x<=T ; x++)	
		cout<< q_star_t[x] << " " ;

	cout<<"\nTotal Iterations : "<<itr<<endl << endl;
	return;
}

//---------------------------------------------
//Before each model run , set all models to 0
void reset_arrays(){
	memset(alpha, 0, sizeof(alpha));
	memset(beta, 0, sizeof(beta));
	memset(gamma, 0, sizeof(gamma));
	memset(delta, 0, sizeof(delta));
	memset(psi, 0, sizeof(psi));
	memset(zeta, 0, sizeof(zeta));
	return;
}

void reset_avgModel(){	//set the model to 0
	memset(A_avg, 0, sizeof(A_avg));
	memset(B_avg, 0, sizeof(B_avg));
	memset(Pi_avg, 0, sizeof(Pi_avg));
	return;
}

//print final averaged models
void diplay_avgmodel() //display the whole model
{
	cout << "State probabilities: " << endl;
	for(int i=1;i<=N;i++)  //display initial state state distributions
	{
        cout << Pi_avg[i] << "\t";
	}
	cout << "\ntransition probabilities: " << endl;

    for(int i=1;i<=N;i++) //display state transition probabilities
	{
        for(int j=1;j<=N;j++)
		{
            cout << A_avg[i][j] << "\t";	
		}
		cout << endl;
	
	}
	cout << "observation symbol probability: " << endl;
	
    for(int i=1;i<=N;i++) //display observation symbol probability distribution
	{
        for(int j=1;j<=M;j++)
		{
            cout << B_avg[i][j] << "\t";
		}
		cout << endl;

	}
	return;
}

//--------------------------------------------

//Storing final models to file
void print_model(ofstream& out) //display the whole model
{
	//cout<< "Printing AVG of all models to output file : "<<endl;
	out<<endl<< "================================================================" <<endl;
	out << "State probabilities: " << endl;
	for(int i=1;i<=N;i++)  //display initial state state distributions
	{
        out << Pi_avg[i] << "\t";
	}
	out << "\ntransition probabilities: " << endl;

    for(int i=1;i<=N;i++) //display state transition probabilities
	{
        for(int j=1;j<=N;j++)
		{
            out << A_avg[i][j] << "\t";	
		}
		out << endl;
	
	}
	out << "observation symbol probability: " << endl;
	
    for(int i=1;i<=N;i++) //display observation symbol probability distribution
	{
        for(int j=1;j<=M;j++)
		{
            out << B_avg[i][j] << "\t";
		}
		out << endl;
	}
	return;
}
void print_avg_a(ofstream& out){	//storing avg models
	 for(int i=1;i<=N;i++) 
	{
        for(int j=1;j<=N;j++)
		{
            out << A_avg[i][j] << "\t";	
		}
		out << endl;
	}
}
void print_avg_b(ofstream& out){
	for(int i=1;i<=N;i++) 
	{
        for(int j=1;j<=M;j++)
		{
            out << B_avg[i][j] << "\t";
		}
		out << endl;
	}
}
void print_avg_pi(ofstream& out){
	for(int i=1;i<=N;i++)  //display initial state state distributions
	{
        out << Pi_avg[i] << "\t";
	}
}



//------------------------------------------------------------------------------------

//Live Testing 
void live_testCi(const char* path){		//Get Ci of live test 
	genrate_Uni(-2 , -1 , path);
	//printtestCi();	
}

void offline_testCi(const char* path){	//get Ci of test data
	genrate_Uni(-2 , 0 , path);
}

void get_test_O(){	//get obs seq for test from codebook and Cis
	double tok_wt[p] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

	for(int i = 0 ; i < T ; i++){			//For each of frame of test.

		double minDist = DBL_MAX ; int minDistIndx = -1 ;

		for(int j = 0; j< K; j++){				//iterate each of CB.
			double temp = 0;

			for(int ci = 0; ci < p ; ci++){			//cur frame
				temp += ( tok_wt[ci] * pow((test_Ci[i][ci] - codebook[j][ci] ),2) ) ;
			}
			//cout<< j+1 << " " << temp << endl;

			if(temp < minDist){
				minDist= temp;
				minDistIndx= j+1;
			}
		}
		O[i+1] = minDistIndx;
		//cout<< i+1 << " " << O[i] << endl;
	}
	return ;

}

//------------------Pipelines------------

void pipeline2obs(){ //full process for getting obs seq from raw samples. 
	get_Uni();	
	dump_Uni();
	Uni_2d();	
	read_cbtxt() ;		//use TA codebook.
	get_obseravtion();	
	dump_obs() ;	
}


void pipeline2obs2(){ 
	get_Uni();	
	dump_Uni();
	Uni_2d();	
	gen_cb(); dump_cb() ; read_cb();	//create own codebook
	get_obseravtion();	
	dump_obs() ;	
}
//-------------------------------------------------------------

int _tmain(int argc, _TCHAR* argv[])
{

	//-------------------------------------------------
	
while(true){
	int choice;
	cout<<"\n------------HMM ASSIGNMENT--------------------\n";
	cout << "\n1 -> Train , 2 -> offline test  , 3 -> Live test : \nEnter a choice :"; cin>> choice;
	while(cin.fail()) {
        cout << "Enter only Interger !!"  << std::endl;
        cin.clear();
        cin.ignore(256,'\n');
        cin >> choice;
     }
	

	//-------------------Training--------------------
	if(choice == 1){
		cout<< " Preprocessing Samples , Creating Codebook from them using LBG, and generating Observation sequence from it..........\n";
		pipeline2obs();
		read_obs() ;
		cout<< "\nObs Seq Generated. Check file Obs Seq.csv , codebook.csv \n";
		cout<<"Now, Training on given observations...\n";

		char filenum[3] ;
		ofstream a_avg ; ofstream b_avg ; ofstream pi_avg ; 
		ofstream out;
	
		int iters = 200;	int utters = 20;

		char output[30] = "output.txt" ;
		out.open(output);

		for(int d = 0 ; d <=9 ; d++){
			char a[30] = "./models/a_avg";			//storing avg of 20 models. 
			char b[30] = "./models/b_avg";   
			char pi[30] = "./models/pi_avg";   
			sprintf(filenum ,"%d", d) ;

			strcat(a, "_") ; strcat(b, "_") ; strcat(pi, "_") ;
			strcat(a, filenum) ; strcat(b, filenum) ;strcat(pi, filenum) ;	
			strcat(a, ".txt") ; strcat(b, ".txt") ;strcat(pi, ".txt") ;
			a_avg.open(a); b_avg.open(b); pi_avg.open(pi);

			reset_avgModel(); //resets avgmodel arrays to 0;
			cout<< "Training Model for digit: " << d << endl;

			for(int u = 1 ; u <= utters ; u++){
				reset_arrays();			//resets all arrays to 0;
				feed_forward_model();

				get_O2(d,u);	       //fills O with current obs seq.				

				cout<< "utterance : " << u << "......" << "  " ;
				RunModel(iters);
				add2avg();	            //sum all models into avg model arrays
			}
			do_avg(utters);		         //get avg of 20 models - now we have A_avg,B_avg,Pi_avg as best model

			cout<<endl;
			out<< endl<< "Avg Model for digit: " << d <<endl ;
			print_model(out);	
			print_avg_a(a_avg);  print_avg_b(b_avg);  print_avg_pi(pi_avg);	//print in respective model files
			a_avg.close() ; b_avg.close() ; pi_avg.close();

			//There wont be any final p_star. Its based per observation.

		}
		cout<< "\n Check file Obs Seq.csv , universe.csv , codebook , avg models in folder models \n";
		out.close();

	}
	//-----------------OFFLINE TESTING-------------------------
	else if(choice == 2){
		read_cbtxt() ;		//read codebook
		int utter ;  ld fp; //forward proc value store
		char filenum[3] ;
		ld C =0 , c = 0; 
		double acc , final_acc ;

		for(int dig = 0 ; dig <= 9 ; dig++){

			for(utter = 1; utter <= 10; utter++){
				//open file
				char path[80] = "./tests/test_";	

				sprintf(filenum ,"%d", dig);    strcat(path, filenum) ;	
				strcat(path, "_") ;
				sprintf(filenum ,"%d", utter);		strcat(path, filenum) ;
				strcat(path, ".txt") ;

				offline_testCi(path);	
				get_test_O() ;	

				cout<< "testing digit :" << dig << " , utterance : "<< utter << endl;
				ld best_p = 0 ; int d = 0;

				for(int m= 0 ; m <= 9 ; m++){ //checking
					get_avg_model(m); 

					cout<<endl<< "Prob of test word being digit "<< m << " = ";

					fp = forward_proc();
					cout << fp << endl;
					if(fp > best_p ){
						best_p = fp;
						d =  m;
					}	
				}
				cout << "Predicted digit : " << d << endl <<endl ;
				if(d == dig){
					c++ ; C++;
				}
			}
			acc = c/10.0;
			cout << "Accuracy for digit " << dig << " = " << acc*100 << "%\n";
			c = 0;
			
		}
		final_acc = C/100.0;
		cout << "\nOverall Accuracy of Model : " << final_acc*100 << "%\n" ;

	}

	//-------------------------------------
	//Live testing
	else if(choice == 3){
		while(true){

			char path[30] = "livetest.txt";	 //generate test
			char cmd[100] = "Recording_Module.exe 3 i.wav ";

			strcat(cmd, path) ;
			system(cmd);

			live_testCi(path);	     //get Ci of test from get_ceps.h
			read_cbtxt() ;		    //read codebook from TA

			get_test_O() ;		   //get obs seq from Ci and codebook
		
			cout << "The observation seq we got for the live test : " << endl;		//print obs seq of test
			printO() ;

			ld fp; ld best_p = 0 ; int dig = 0;

			for(int m= 0 ; m <= 9 ; m++){		//digit
				get_avg_model(m);
				cout<< endl<< "Probability of spoken word being the digit : "<< m << " = "<< endl ;
				fp = forward_proc();
				cout << fp << endl;

				if(fp > best_p ){
					best_p = fp;
					dig =  m;
				}
			}
			cout << "The Predicted digit : " << dig << endl <<endl ;

			cout<< "Press 1 to run test again else 2 to exit. ";
			int ch; cin >> ch;
			while(cin.fail()) {
				cout << "Enter only Interger !!"  << std::endl;
				cin.clear();
				cin.ignore(256,'\n');
				cin >> ch;
			 }
			if(ch != 1)break;
		}
	}

	//------------------------------
	else{
		cout << "Wrong choice !! Exiting !!!";
		return 0;
	}

	}

	//-----------------------------------------------------
	
}

