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
#include "Controller.h"
#include "Cooling.h"
#include "Zone.h"

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
	ParallelComponentsComponentsController* control = new ParallelComponentsComponentsController();
	ParallelComponentsComponentsCoolingSystem* cooling =new ParallelComponentsComponentsCoolingSystem();
	ParallelComponentsComponentsZone* zone = new ParallelComponentsComponentsZone();
	fmus.push_back(control);
	fmus.push_back(cooling);
	fmus.push_back(zone);
	// set the start and end time of simulation
	double Tstart = 0;
	double Tend = 1000;
	// fixed step size
	double dt = 0.1; 
	double Tnext = Tend;
	double time  = Tstart;
	// Mapping variables 
	double val;
	// - USER DEFINES - END!

	//Initialize all fmus
	int size = fmus.size();
	#pragma omp parallel for
	for(int i = 0;i<size;i++){
		
		FMI* fmi = fmus[i];
		// Set initial value for time
		fmi->_fmi2SetTime(fmi->c,time);
		// Setup
		fmi->_fmi2SetupExperiment(fmi->c,fmi2False,0.0,Tstart,fmi2True,Tend);
		fmi->q = new fmi2Real[fmi->numStates];
		fmi->q_pre = new fmi2Real[fmi->numStates];
		fmi->q_nominal = new fmi2Real[fmi->numStates];
		fmi->dq = new fmi2Real[fmi->numStates];
		fmi->z = new fmi2Real[fmi->numEvents];
		fmi->pre = new fmi2Real[fmi->numEvents];
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
		// Map inputs and outputs here ########################################################## - USER DEFINES - START!
		// Temperature sensor to control
		val = zone->get_y();
		control->set_u(val);
		// Control sets air_flow of cooling
		val = control->get_y();		
		cooling->set_u(val);
		// Fan blows air in the thermal zone
		// Cooling to Zone
		val = cooling->get_outlet_m_flow();
		zone->set_inlet1_m_flow(val);
		val = cooling->get_outlet_forward_T() ;
		zone->set_inlet1_forward_T(val);
		val = cooling->get_outlet_forward_X_w() ;
		zone->set_inlet1_forward_X_w(val);
		// Zone to Cooling
		val = zone->get_outlet1_m_flow();
		cooling->set_inlet1_m_flow(val);
		val = zone->get_outlet1_forward_T();
		cooling->set_inlet1_forward_T(val);
		val = zone->get_outlet1_forward_X_w();
		cooling->set_inlet1_forward_X_w(val);
		// Calculate next states and finalize the step ########################################## - USER DEFINES - END!
		// Rolling back if there are one or more events that overshoots the tolerence while Z is changing sign 
		if (rollback == true){
			time=time-h;
			double t_next = Tend;
			//#pragma omp parallel for
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
				#pragma omp parallel for
				for(int j =0;j<fmi->numEvents;j++){
					fmi->z[j]=fmi->pre[j];	
				}
				#pragma omp parallel for
				for(int j =0;j<fmi->numStates;j++){
					fmi->q[j]=fmi->q_pre[j];	
				}
				fmi->_fmi2SetContinuousStates(fmi->c,fmi->q_pre,fmi->numStates);
			}
		}
		// If no rollback then execute events for the Z values that are lower than tolerence
		if (rollback == false){
			#pragma omp parallel for
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
	//cout<<"Time:	"<<time<<"	Room Temp:	"<<zone->get_senTemRoo_T()<<"	Y:	"<<zone->get_y()<<"	U	"<<control->get_u()<<endl;
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

