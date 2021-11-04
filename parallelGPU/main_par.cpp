#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include "omp.h"
#include "time.h"
#include <cstdlib>
#include <stdlib.h>
#include <dlfcn.h>
#include "assert.h"
#include "omp.h"
#include <vector>

using namespace std;
// Tank class	
class Tank {
public:
	double x;
  	double a;
  	double u;
	Tank(){
		x = 10;
		a = 0.2;
	}
	double derivative(double u, double x, double a){
	
		double derx=u-a*x;
		return derx;
	}
	~Tank(){}
};
	
int main(int argc, char *argv[]) {
	
	// Read the arguments
	// The number of Tanks in each FMU    
    	std::istringstream is( argv[1] );
        int numOfTanks;
        if (is >> numOfTanks)
        {
            // Conversion successful
        } 
	// The number of Replications
	std::istringstream iss( argv[2] );
        int reps;
        if (iss >> reps)
        {
            // Conversion successful
        }
    	// The number of Cores
	std::istringstream issss( argv[3] );
        int numOfCPUs;
        if (issss >> numOfCPUs)
        {
            // Conversion successful
        }  
   	 // Set the number of threads
	omp_set_num_threads(numOfCPUs);      
   	// For each replication
	for (int n = 0;n<reps;n++){
		
		// Record execution time	
		double start = omp_get_wtime();
		// Create FMU set
		vector<vector<Tank*>> fmus;
		//#pragma omp parallel for
		for(int i = 0;i<numOfCPUs;i++){
			// Create Tanks vector
			vector<Tank*> tanks;
			for(int j = 0;j<numOfTanks;j++){
	
				Tank *t = new Tank();
				tanks.push_back(t);		
			}
			fmus.push_back(tanks);
		}
		// set the start and end time of simulation
		double Tstart = 0;
		double Tend = 10000.0;
		// fixed step size of 10 milli-seconds
		double dt = 0.1; 
		double time  = Tstart;
		// Forward Euler
	   	while (time < Tend)
	   	{	
	   		// advance time
			time = time + dt;
			// Map inputs and outputs of Tanks
			#pragma omp parallel for
			for(int z = 0;z<numOfCPUs;z++){
				
				// Here we will add GPU parallelization - #pragma acc -directive suc as parallel XXXXXXXXXXXXXXXXXXXXXXXXXXXX 
		   		for(int i = 0;i<numOfTanks;i++){
				
					// Not the First Tank
					if (i > 0){
						
						// Get input values from the previous Tank
						double val = fmus[z][i-1]->x*fmus[z][i-1]->a;
						fmus[z][i]->u=val;		
					}
					else {
						// First batch
						if (z == 0){
							
							fmus[z][i]->u=1.0;
						}
						else{
							double val = fmus[z-1][numOfTanks-1]->x*fmus[z-1][numOfTanks-1]->a;
							fmus[z][i]->u=val;
						}
					}
				}
			}
			#pragma omp parallel for
			for(int z = 0;z<numOfCPUs;z++){
				// Here we will add GPU parallelization - #pragma acc -directive suc as parallel XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
				for(int i = 0;i<numOfTanks;i++){

					// compute derivatives
					double derx = fmus[z][i]->derivative(fmus[z][i]->u,fmus[z][i]->x,fmus[z][i]->a);
					fmus[z][i]->x = fmus[z][i]->x + dt*derx;
				}
			}
			// print results
			//cout<<" "<<time<<" Tank1"<<" "<<fmus[0][0]->x<<" Tank2"<<" "<<fmus[0][1]->x<<endl;
		}
		
		// cleanup
		for(int i = 0;i<numOfCPUs;i++){
			
			for(int j = 0;j<numOfTanks;j++){
				
				delete fmus[i][j];
			}
		}

		double end = omp_get_wtime();
		double seconds = end-start;	 
		cout <<seconds<< endl;
	}	
}

