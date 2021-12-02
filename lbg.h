
void KmeansClustering(int size2, int k ){					//Will be used For k-means as well for LBG method. 
	int m = 0;           // No of iterations
	int i = 0 , j = 0;
	double old_dist = 0 , cur_dist = 1;     // holds current and prev distortion
	double tok_wt[p] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0}; //tokuhara weights
	double dlta = 0.00001 ;       // exit condition 
	int vec =0, cbvec = 0;		//index to iterate in Universe and codebook respectively.
	double dist,tok_dist = 0;	// to caculate tokuhara distance

	double min_dist = DBL_MAX ; // holds nearest neighbour distance 

	//step 2-6
	int ch = 0;

	while( abs(old_dist - cur_dist) > dlta ){	  //exit condition
		old_dist = cur_dist ;					// make prev iteration's distance as old distance
	
		vector<int> region_array(size2);
		int vec_count[K] = {0};			    //stores count of vectors in each region.
	

		for(vec = 0; vec < size2 ; vec++){		//for each vector from universe , calculate distance from each of codebook vectors
			min_dist = DBL_MAX;
			int region = 0 ;				   // holds region to assign for a vector
			for(cbvec = 0 ; cbvec <k ; cbvec++){ 
				tok_dist = 0;	
				for(j = 0 ; j < p ; j++){
					dist = Uni2d[vec][j] - codebook[cbvec][j] ;
					tok_dist += tok_wt[j]* (dist*dist);
				}			
				if(tok_dist < min_dist){		// if current vector has least distance with current codebook vector,	
					region = cbvec ;			// set region as index of current codebook vector (cbvec),
					min_dist = tok_dist;		// and set that distance as least.
				}	
			}
			region_array[vec] = region;  // after min dist from codebook is found, assign training vector to that region

			vec_count[region]++ ;		// and increase count of vectors in that region
			
		}
		//step 3 and 4
		double mean_points[K][p] = {{0}} ;    // will store new centroids for each region
		cur_dist = 0; dist = 0; tok_dist = 0;

		for(cbvec = 0 ; cbvec < k ; cbvec++){    // for every region
			for(vec = 0 ; vec < size2 ; vec++){		//check which vector from universe falls in cur region
				if(region_array[vec] != cbvec)continue;	//if it doesnt belong to cur region, move to next vector

				//step 3  - Getting distortion from all regions and averaging by dividing by size2 
				for(j = 0 ; j < p ; j++){	   // calculate its distortion .
					dist = codebook[cbvec][j] - Uni2d[vec][j] ;
					tok_dist += tok_wt[j] * (dist*dist);  //summing tokuhara distance for all points , for all vectors

					//step 4 - Summing all vector points of a  region to average later and get centroid
					mean_points[cbvec][j] += Uni2d[vec][j];
				}
			}	 
		}
		//tok_dist now contains distortion for all size2 vectors 
		cur_dist = tok_dist/size2 ; // avg tokuhara dist and store in cur_dist - step 3 done for current iteration

		//For updating centroids - dividing mean_points sum of all vector in each region
		//by no of vectors in that region to get updated codebook/centroid
		for(cbvec = 0 ; cbvec < k ; cbvec++){
			for(j = 0 ; j < p ; j++){
				//cout<< vec_count[cbvec] << " ";
				codebook[cbvec][j] = mean_points[cbvec][j]/vec_count[cbvec];
			}
		}
		//Now we have updated codebook - step 4 done for current iteration
		
		m++ ;		//step 5
	} //end of while loop
}


void get_codebook()	//No use of , d , u.
{
	
	int size2 = unisize;             // size2 = size of universe	
	//cout << size2;
	//------------------------------------------------------------------------//

	int k = 1;          //k of k-means for each iteration of LBG
	int vec = 0, j=0;  //vec to loop universe and j to loop all points of vector

	//step 1 - initialize codeboook 
	for(vec = 0; vec < size2; vec++){    //iterate all vectors
		for(j=0; j < p; j++){			// iterate all points of a vector
			codebook[0][j] += Uni2d[vec][j];		 //1st vector of codebook will be our initialised centroid vector
		}
	}
	for(j=0; j < p; j++){			
		codebook[0][j] = codebook[0][j]/size2 ;   //get average which will be our initialised centroid vector 
	}

	double epsilon = 0.03 ; 
	while(k < K){    //iterations runs till k=8
		//step 2 - Splitting vectors- We split from last centroid vector(k-1)th to 1st centroid vector(0)th 
		// Split them in 2 vectors, so, for k=1 , we get 2 vectors ,and 4 in next iteartion and then 8 for next iteration.
		//and store the 2 new vectors in index "2*vec+1" and "2*vec" for vector of index "vec".
		
		for(vec = k-1 ; vec >= 0 ; vec--){  // iterating codebook vectors which are centroid now, from last centroid to first
			for(j = 0 ; j < p ; j++){
				codebook[2*vec+1][j] = codebook[vec][j]*(1+epsilon); // modifying and storing at respective index as explained above.
				codebook[2*vec][j] = codebook[vec][j]*(1-epsilon);		
				//for k=1 , we now have 2 vectors at index 0 and 1. Index 1 is created first ,
				//then index 0 , so old vector is updated at index 0.(no need of old vector)
				//So for each iteration , codebook is doubled.
			}
		}
		k = k*2;	//telling k means about doubled size2 of codebook. 

		//Now , calling k means on this codebook for given k. 
		KmeansClustering(size2,k) ;    //k = new size2 of codebook AKA no of rows AKA no of clusters.
	}
	return ;
}

