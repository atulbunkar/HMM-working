

void printcb(){
	for(int i=0; i < 32; i++){
		cout << i << endl;
		for(int j=0; j < 12; j++){
			cout<< codebook[i][j] << " " ;
		}
		cout<< endl;
	}
}

void printO(){
	for(int i =1; i<= T; i++)
		cout<< O[i] << " " ;
	cout<<endl << endl;
	return;
}

void printtestCi(){
	for(int frame = 0; frame < T ; frame++){
		for(int x=0; x < p ; x++){	
			cout << test_Ci[frame][x] << "  "; 
		}
		cout << endl;
	}
	return;


}

void get_avg_a(int d){	
	char num[3];
	char fullname[100] = "./models/a_avg_";
	ld number; 
	int row=1; int idx = 1;

	sprintf(num ,"%d", d);
	strcat(fullname, num);
	strcat(fullname,".txt");
	//printf("\n %s", fullname );

    ifstream file(fullname);
	if(!file)cout<< "\n cant open";
	string input;

	while(getline(file, input)){
		istringstream iss( input );
		idx = 1; 
		while( iss >> number ){
			A[row][idx++] = number;
		}
		row++;
	}
	file.close();
	return  ;

}

void get_avg_b( int d){	
	char num[3];
	ld number;
	char fullname[100] = "./models/b_avg_";

	sprintf(num ,"%d", d);
	strcat(fullname, num);
	strcat(fullname,".txt");
	//printf("\n %s", fullname );

    ifstream file(fullname);
	if(!file)cout<< "\n cant open";
	string input; int idx = 1; int row = 1;

	while(getline(file, input)){
		istringstream iss( input );
		idx = 1; 
		while( iss >> number ){
			B[row][idx++] = number;
		}
		row++;
	}
	file.close();
	return ;
}

void get_avg_pi(int d){
	char num[3];
	char fullname[100] = "./models/pi_avg_";
	 int number;

	sprintf(num ,"%d", d);
	strcat(fullname, num);
	strcat(fullname,".txt");
	//printf("\n %s", fullname );

    ifstream file(fullname);
	if(!file)cout<< "\n cant open";
	string input;

	getline(file, input);
	istringstream iss( input );
	int idx = 1;
	while( iss >> number ){
		Pi[idx++] = number;
	}
	file.close();
	return ;
}

void get_avg_model( int d){
	get_avg_pi(d) ;
	get_avg_a(d) ;
	get_avg_b(d) ;
}



