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
#include "ParallelComponentsComponentsOneZoneOneCooler1.h"
#include "ParallelComponentsComponentsOneZoneOneCooler2.h"
#include "ParallelComponentsComponentsOneZoneOneCooler3.h"
#include "ParallelComponentsComponentsOneZoneOneCooler4.h"

using namespace std;
	
int main(int argc, char *argv[]) {
	
	////// Read the arguments
	// The number of CPUs
	std::istringstream iss( argv[1] );
        int numOfCPUs;
        if (iss >> numOfCPUs){
            // Conversion successful
        } 
    	// Set the number of threads
	omp_set_num_threads(numOfCPUs);      
		
	// Record execution time	
	double start = omp_get_wtime();
	// Create FMI object
	vector<FMI*> fmus;
	// Create FMI objects list - USER DEFINES - START!
	ParallelComponentsComponentsOneZoneOneCooler1* room1 = new ParallelComponentsComponentsOneZoneOneCooler1();
	ParallelComponentsComponentsOneZoneOneCooler2* room2 = new ParallelComponentsComponentsOneZoneOneCooler2();
	ParallelComponentsComponentsOneZoneOneCooler3* room3 = new ParallelComponentsComponentsOneZoneOneCooler3();
	ParallelComponentsComponentsOneZoneOneCooler4* room4 = new ParallelComponentsComponentsOneZoneOneCooler4();
	fmus.push_back(room1);
	fmus.push_back(room2);
	fmus.push_back(room3);
	fmus.push_back(room4);
	// set the start and end time of simulation
	double Tstart = 0;
	double Tend = 10000;
	// fixed step size
	double dt = 0.1; 
	double Tnext = Tend;
	double time  = Tstart;
	// - USER DEFINES - END!
	// Mapping variables 
	double val;

	//Initialize all fmus
	int size = fmus.size();
	//#pragma omp parallel for
	for(int i = 0;i<size;i++){
		
		FMI* fmi = fmus[i];
		fmi->q = new fmi2Real[fmi->numStates];
		fmi->q_pre = new fmi2Real[fmi->numStates];
		fmi->q_nominal = new fmi2Real[fmi->numStates];
		fmi->dq = new fmi2Real[fmi->numStates];
		fmi->z = new fmi2Real[fmi->numEvents];
		fmi->pre = new fmi2Real[fmi->numEvents];
		// Set initial value for time
		fmi->_fmi2SetTime(fmi->c,time);
		// Setup
		fmi->_fmi2SetupExperiment(fmi->c,fmi2False,0.0,Tstart,fmi2True,Tend);
		// Initiate all fmus
		fmi->init(fmi->q, time);
	}
   	while (time < Tend){	
		double h = 0;
		#pragma omp parallel for
   		for(int i = 0;i<size;i++){
			FMI* fmi = fmus[i];
			// Compute derivatives
			fmi->der_func(fmi->q,fmi->dq);
		}
		// Time advance ##########################################################
		h = std::min(dt, Tnext);
		time = time+h;

		// Map inputs and outputs here ########################################################## - USER DEFINES - START!
		// T 1-2
		val = room1->get_y();
		room2->set_u(val);
		val = room2->get_y();
		room1->set_u(val);
		// T 1-3
		val = room1->get_y();
		room3->set_u1(val);
		val = room3->get_y();
		room1->set_u1(val);
		// T 2-4
		val = room2->get_y();
		room4->set_u1(val);
		val = room4->get_y();
		room2->set_u1(val);
		// T 3-4
		val = room3->get_y();
		room4->set_u(val);
		val = room4->get_y();
		room3->set_u(val);
		// Calculate next states and finalize the step ############################################ - USER DEFINES - END!
		#pragma omp parallel for				
		for(int i = 0;i<size;i++){
			FMI* fmi = fmus[i];
			// Set time
			fmi->_fmi2SetTime(fmi->c,time);			
			for(int j = 0; j<fmi->numStates;j++){
				
				fmi->q_pre[j]=fmi->q[j];	
				fmi->q[j]=fmi->q[j]+h*fmi->dq[j];
			}
			fmi->_fmi2SetContinuousStates(fmi->c,fmi->q,fmi->numStates);
			fmi->postStep(fmi->q,fmi->z,fmi->pre,Tend, time, h);
		}
		// Check if there is an event that has a zero crossing function value that overshoots the tolerence
		Tnext = Tend;
		bool rollback = false;
		for(int i = 0;i<size;i++){
			FMI* fmi = fmus[i];
			if (fmi->event == true && fmi->cont_time_mode == true){ 
				rollback = true;
			}
		}

		// Rolling back if there are one or more events that overshoots the tolerence while Z is changing sign 
		if (rollback == true){
			time=time-h;
			double t_next = Tend;
			#pragma omp parallel for
			for(int i = 0;i<size;i++){
				FMI* fmi = fmus[i];
				// Put all in continous state mode
				fmi->cont_time_mode = true;
				if (fmi->next_event_time != 0 && fmi->next_event_time < t_next){
						
					t_next = fmi->next_event_time;
					// Assign new event time
					Tnext=t_next;
				}
				// Rollback state variables, Z, and time
				fmi->event = false;
				fmi->_fmi2SetTime(fmi->c,time);
				for(int j =0;j<fmi->numEvents;j++){
					fmi->z[j]=fmi->pre[j];	
				}
				for(int j =0;j<fmi->numStates;j++){
					fmi->q[j]=fmi->q_pre[j];	
				}
				fmi->_fmi2SetContinuousStates(fmi->c,fmi->q_pre,fmi->numStates);
			}
		}
		// If no rollback then execute events for the zero-crossing values that are lower than tolerence
		if (rollback == false){
			for(int i = 0;i<size;i++){
				FMI* fmi = fmus[i];
				// Do events
				if (fmi->cont_time_mode == false){
					fmi->doEvents(fmi->q,fmi->q_nominal, Tend);
					Tnext = Tend;
				}
			}
		}							
	// Output
	//cout<<"Time:	"<<time<<"	Room1 Temp:	"<<room1->get_senTemRoo_T()<<"	Room2 Temp:	"<<room2->get_senTemRoo_T()<<"	Room3 Temp:	"<<room3->get_senTemRoo_T()<<"	Room4 Temp:	"<<room2->get_senTemRoo_T()<<endl;
	}
	// cleanup
	for(int i = 0;i<size;i++){
		
		delete fmus[i];
	
	}
	// Measure execution time
	double end = omp_get_wtime();
	double seconds = end-start;
	cout <<seconds<< endl;
}

