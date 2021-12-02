int skip_line(FILE *fp) //  function to skip lines, i.e move cursor to next line !
{
    int c;

    while (c = fgetc(fp), c != '\n' && c != EOF);

    return c;
} 

void get_samples(FILE* fp , double *samples ){  // Read 320 samples from test.txt for step 1
	char line_buffer[20];
	double Xi; int c = 0;
	double sum =0.0; 
	int j = 0 , k =0;
	
	while(!feof(fp)){          
		fgets(line_buffer, 20, fp);
		Xi = atof(line_buffer);
		if(Xi != 0.0){
			samples[c] = Xi;
			c++;
		}
	}
}

void calc_ri(double *ri , double *samples){ // calculates autocorrelation of data and stores in ri
	int j = 0 , k =0;
	double sum = 0.0;
	for(k; k <= p ; k++){
		sum = 0.0;
		for(j=0 ; j < F-k ; j++){
			sum += samples[j]*samples[j+k]*1.0; 
		}
		ri[k] = sum;
	}

}

void calc_ai(double *ri , double *ai){ // calculates ai's using durbin algo and stores in ai
	double e[p+1] , a[p+1][p+1] ,k[p+1] ;
	int i= 0 ,j = 0 ; double sum = 0.0;
	e[0] = ri[0];
	for(i=1 ; i <= p ; i++){
		sum = 0;
		for(j = 1 ; j <= i-1 ; j++){
			sum+= a[i-1][j]*ri[i-j]*1.0;
		}
		k[i] = (ri[i]-sum)*1.0/e[i-1];
		a[i][i] = k[i];
		for(j =1 ; j <= i-1 ;j++){
			a[i][j] = a[i-1][j] - k[i]*a[i-1][i-j]*1.0;
		}
		e[i] =(1-(k[i]*k[i])) * e[i-1] *1.0;
	}
	for(i=1;  i <= p ; i++){
		ai[i] = a[p][i];
	}

}

void calc_ci(double *ri ,double *ai , double *ci){  // calculates cepstral's coefficients  ci's and stores in ci
	double sigma = ri[0]; int m ,k ; double sum = 0; 
	
	ci[0] = log(sigma*sigma) ;
	for(m = 1; m <= p ; m++){
		sum = 0;
		for(k = 1; k <= m-1 ; k++){
			sum += (k*ci[k]*ai[m-k] )*1.0/m;
		}
		ci[m] = ai[m] + sum ; 
	}
}
 
void print(double *arr , int start ,int end){		// print array elements
	int x ;
	for(x =start; x <=end; x++)
		printf("%lf \t" , arr[x]);
	printf("\n");
}

void raised_sine(double *ci){ //applying the raised sine window on ci
	int m;
	double pi = 3.141 ;

	for(m = 1 ; m <= p ; m++){
		ci[m] = ci[m] * (1.0 +  (p/2.0)* sin( pi*m*1.0/p )) ;   
	}
	
}

int preprocess(FILE* fp_raw , int u){
	double dc_bias = 0; int sampleCount = 0;
	int frame_size = 100 ; double norm_mult;
	double Xi=0 ;
	double peak_Xi=0 ; int line_skips;
	char line_buffer[20] ; int new_peak = 5000 ;         // new_peak will be the peak signal after normalization.
	int i =0 ; double xi = 0.0, norm_xi = 0.0;
	int fskip = 60;										//while live testing to get rid of massive initial DC shift , wait for 1 sec before speaking word.

	if(u==-1){
		line_skips = fskip*F ;
		while(line_skips--){
			skip_line(fp_raw) ; 
		}
	}

	while(!feof(fp_raw)){
		fgets(line_buffer, 20, fp_raw);		//getting dc bias
			double x = atof(line_buffer);
			dc_bias += x;
			sampleCount++;
	}
	dc_bias /= sampleCount;
	//cout<< "dc " << dc_bias;
	
	fseek(fp_raw , 0 , SEEK_SET);
	if(u==-1){
		line_skips =fskip*F ;
		while(line_skips--){
			skip_line(fp_raw) ; 
		}
	}

	// Normalization                   
	while(!feof(fp_raw)){ 
		fgets(line_buffer, 20, fp_raw);						 // Here, I found the max amplitude signal from the data.
			Xi = atof(line_buffer);
			if(abs(Xi) > peak_Xi ){
				peak_Xi = abs(Xi);
			}    
	 }

	norm_mult = new_peak*1.0/abs(peak_Xi) ;  // This gives us the normalization factor. 

	fseek(fp_raw , 0 , SEEK_SET);				
	if(u==-1){
		line_skips = fskip*F ;
		while(line_skips--){
			skip_line(fp_raw) ; 
		}
	}
	
	while(!feof(fp_raw)){                        // I am transforming the raw data here by subtracting from each point the dc_bias , then multiplying them with 
		fgets(line_buffer, 20, fp_raw);			 // the normalization factor.
		xi = atof(line_buffer);
		norm_xi = ((xi-dc_bias) * norm_mult*1.0 );      
		if(norm_xi != 0)                                
			norm_data[i++] = norm_xi;        // Storing the normalised data in array "norm_data".
			//cout << norm_xi << " ";
		}
	return i;								// returns the number of data points present in normalized data !
}

int find_start(int s){					//finds stable frame  s->size of preprop data.
	int start = 0;
	int st = 0;      
	double cur_Xi = 0, energy = 0;
	double max_energy = 0 ;

	while(start < s){                                  // This loops in normalised data.
		cur_Xi =  norm_data[start++];

		energy += cur_Xi*cur_Xi;                     // calculating energy for current frame
		
		if(start%F == 0){						   // check end of frame.
			energy = energy*1.0/F ;			      // dividing by num of data .  														      
			if(energy > max_energy){			 // found max energy frame
				st = start-F ;					// start is at end of frame so st has start of frame sample point's index.
				max_energy = energy;		
			}
		}
	}
	int x = int(T/3);	
	st = int(st - x*slide);		//taking T/2 frame behind max energy frame as starting frame.
	return st ;	
}

ld ambient_energy_calc(){
	ld cur_Xi=0.0 , prev_Xi=0.0;                 // cur_Xi is the current data point , prev_Xi is prev data point.
	int sample_num =0; 
	ld ambient_energy = 0 ,avg_ambient_energy=0 , ambient_zcr =0;
	int frame_size = 320;

	while(sample_num<=3*F){               // This loop iterates for 3 frames of the start of the word which is silence and gets its energy and zcr.
		cur_Xi = norm_data[sample_num];
		if(sample_num>0){
			ambient_energy += cur_Xi*cur_Xi;	                    // storing sum of square of Xi , will get average after frame ends. 
		}
		prev_Xi = cur_Xi;										// this sets cur Xi as prev for next iteration.
		if(sample_num%F ==0){                                 // checks end of frame , each time data crosses 100 points 
			avg_ambient_energy += ambient_energy/F;     // Summing the avg ambient energy of all frames , later will divide by num of frames i.e 10. so we will get energy of silence .
			ambient_energy = 0;	                                 // After every frame ends , we set energy back to 0 for next iteration.
		}
		sample_num++;
	}
	avg_ambient_energy = avg_ambient_energy/3 ;           // We summed avg energy for 10 frames so dividing to get their avg energy.
	
	return avg_ambient_energy;
}


int find_start2(int s){					//finds stable frame  s->size of preprop data.
	int start = 0;
	int st = 0;      
	double cur_Xi = 0, energy = 1;
	double prev_energy = 10000000000 ;
	ld avg_ambient_energy = ambient_energy_calc();

	while(start < s){                                  // This loops in normalised data.
		cur_Xi =  norm_data[start++];
		energy += cur_Xi*cur_Xi;                     // calculating energy for current frame	
		if(start%F == 0){						   // check end of frame.
			energy = energy*1.0/F ;			      // dividing by num of data .  														      
			if(energy > 2*avg_ambient_energy){			 // found max energy frame
				st = start;
				return st;							//returns at first find.
			}
			prev_energy = energy;
		}
	}
	return st ;	
}

void genrate_Uni(int d , int u ,const char* path){
	//this gives Ci's for LBG.
	double ri[p+1] ;					
	double ai[p+1] ;				  // for storing single frame ai,ci,ri's
	double ci[p+1];
	double samples[F+1];			//stores 320 samples , i.e one frame.
	
	int  j,frame , x =0 ;

	FILE* fpp;
	fpp = fopen(path , "r");
	if (fpp == NULL){						 // checking if file has opened successfully r not.
		printf("error opening file - \n");
	}
	int s = preprocess(fpp,u);	// preprocesses the file and stores in norm_data , returns size(index) till which values are stored.
	fclose(fpp);

	int start = find_start(s);			// finds the start datapoint of stable word
	if(start<0)start = 0;				

	for(frame = 0; frame < T ; frame++){	// runs for T frames 

		for(j =0 ; j < F; j++){				//gives one frame of size F=320.
			samples[j] = norm_data[start+j];	
		}
		start = start+slide;    // now start is moving by 80.- FIX TO SLIDING WINDOW

		calc_ri(ri,samples);		
		
		calc_ai(ri,ai);				
		
		calc_ci(ri,ai,ci);			
		
		raised_sine(ci);			
	
		if(d == -2 ){					//while testing
			for(x=0; x < p ; x++){	
				test_Ci[frame][x] = ci[x+1]; 
			}
			continue;
		}

		for(x=1; x <= p ; x++){			//while training
			Uni[d][u][frame][x] = ci[x];	
			check++;
		}
	}
	return;
}



