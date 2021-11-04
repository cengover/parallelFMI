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
#include "fmi2Functions.h"
#include "fmi2FunctionTypes.h"
#include "fmi2TypesPlatform.h"
#include "assert.h"
#include "omp.h"
#include <vector>

using namespace std;
	
	class FMI {

public:
		/**
		 * This constructs a wrapper around an FMI. The constructor
		 * must be provided with the FMI's GUID, the number of state variables,
		 * number of event indicators, and the path to the .so file
		 * that contains the FMI functions for this model.
		 */	
		 // Create the FMI component
		fmi2CallbackFunctions* callbackFuncs;
		fmi2Component m;
		fmi2Real* fmi_valu;
		fmi2Real* fmi_valx;
		fmi2Real* fmi_vala;
		fmi2Real* fmi_valder;
		
		FMI(const char* modelname,
			const char* guid,
			int num_state_variables,
			int num_event_indicators,
			const char* shared_lib_name,
			const double tolerance = 1E-8);
		// Pointer to the FMI interface
		fmi2Component (*_fmi2Instantiate)(fmi2String, fmi2Type,
				fmi2String, fmi2String, const fmi2CallbackFunctions*,
				fmi2Boolean, fmi2Boolean);
		void (*_fmi2FreeInstance)(fmi2Component);
		fmi2Status (*_fmi2SetupExperiment)(fmi2Component, fmi2Boolean,
				fmi2Real, fmi2Real, fmi2Boolean, fmi2Real);
		fmi2Status (*_fmi2EnterInitializationMode)(fmi2Component);
		fmi2Status (*_fmi2ExitInitializationMode)(fmi2Component);
		fmi2Status (*_fmi2GetReal)(fmi2Component, const fmi2ValueReference*, size_t, fmi2Real*);
		fmi2Status (*_fmi2GetInteger)(fmi2Component, const fmi2ValueReference*, size_t, fmi2Integer*);
		fmi2Status (*_fmi2GetBoolean)(fmi2Component, const fmi2ValueReference*, size_t, fmi2Boolean*);
		fmi2Status (*_fmi2GetString)(fmi2Component, const fmi2ValueReference*, size_t, fmi2String*);
		fmi2Status (*_fmi2SetReal)(fmi2Component, const fmi2ValueReference*, size_t, const fmi2Real*);
		fmi2Status (*_fmi2SetInteger)(fmi2Component, const fmi2ValueReference*, size_t, const fmi2Integer*);
		fmi2Status (*_fmi2SetBoolean)(fmi2Component, const fmi2ValueReference*, size_t, const fmi2Boolean*);
		fmi2Status (*_fmi2SetString)(fmi2Component, const fmi2ValueReference*, size_t, const fmi2String*);
		fmi2Status (*_fmi2EnterEventMode)(fmi2Component);
		fmi2Status (*_fmi2NewDiscreteStates)(fmi2Component,fmi2EventInfo*);
		fmi2Status (*_fmi2EnterContinuousTimeMode)(fmi2Component);
		fmi2Status (*_fmi2CompletedIntegratorStep)(fmi2Component, fmi2Boolean, fmi2Boolean*, fmi2Boolean*);
		fmi2Status (*_fmi2SetTime)(fmi2Component, fmi2Real);
		fmi2Status (*_fmi2SetContinuousStates)(fmi2Component, const fmi2Real*, size_t);
		fmi2Status (*_fmi2GetDerivatives)(fmi2Component, fmi2Real*, size_t);
		fmi2Status (*_fmi2GetEventIndicators)(fmi2Component, fmi2Real*, size_t);
		fmi2Status (*_fmi2GetContinuousStates)(fmi2Component, fmi2Real*, size_t);
		// so library handle
		void* so_hndl;

		static void fmilogger(
			fmi2ComponentEnvironment componentEnvironment,
			fmi2String instanceName,
			fmi2Status status,
			fmi2String category,
			fmi2String message, ...)
		{
			std::cerr << message << std::endl;
		}
	~FMI()
	{
		_fmi2FreeInstance(m);
		delete callbackFuncs;
		dlclose(so_hndl);
	}

};

FMI::FMI(const char* modelname,
			const char* guid,
			int num_state_variables,
			int num_event_indicators,
			const char* so_file_name,
			const double tolerance)
{
	so_hndl = dlopen(so_file_name, RTLD_LAZY);
	if (!so_hndl)
	{
		cout<<"Could not load so file"<< endl;
    }
	
	*(void**) (&_fmi2Instantiate) = dlsym(so_hndl,"fmi2Instantiate");
	*(void**) (&_fmi2FreeInstance) = dlsym(so_hndl,"fmi2FreeInstance");
	*(void**) (&_fmi2SetupExperiment) = dlsym(so_hndl,"fmi2SetupExperiment");
	*(void**) (&_fmi2EnterInitializationMode) = dlsym(so_hndl,"fmi2EnterInitializationMode");
	*(void**) (&_fmi2ExitInitializationMode) = dlsym(so_hndl,"fmi2ExitInitializationMode");
	*(void**) (&_fmi2GetReal) = dlsym(so_hndl,"fmi2GetReal");
	*(void**) (&_fmi2GetInteger) = dlsym(so_hndl,"fmi2GetInteger");
	*(void**) (&_fmi2GetBoolean) = dlsym(so_hndl,"fmi2GetBoolean");
	*(void**) (&_fmi2GetString) = dlsym(so_hndl,"fmi2GetString");
	*(void**) (&_fmi2SetReal) = dlsym(so_hndl,"fmi2SetReal");
	*(void**) (&_fmi2SetInteger) = dlsym(so_hndl,"fmi2SetInteger");
	*(void**) (&_fmi2SetBoolean) = dlsym(so_hndl,"fmi2SetBoolean");
	*(void**) (&_fmi2SetString) = dlsym(so_hndl,"fmi2SetString");
	*(void**) (&_fmi2EnterEventMode) = dlsym(so_hndl,"fmi2EnterEventMode");
	*(void**) (&_fmi2NewDiscreteStates) = dlsym(so_hndl,"fmi2NewDiscreteStates");
	*(void**) (&_fmi2EnterContinuousTimeMode) = dlsym(so_hndl,"fmi2EnterContinuousTimeMode");	
	*(void**) (&_fmi2CompletedIntegratorStep) = dlsym(so_hndl,"fmi2CompletedIntegratorStep");
	*(void**) (&_fmi2SetTime) = dlsym(so_hndl,"fmi2SetTime");
	*(void**) (&_fmi2SetContinuousStates) = dlsym(so_hndl,"fmi2SetContinuousStates");
	*(void**) (&_fmi2GetDerivatives) = dlsym(so_hndl,"fmi2GetDerivatives");
	*(void**) (&_fmi2GetEventIndicators) = dlsym(so_hndl,"fmi2GetEventIndicators");
	*(void**) (&_fmi2GetContinuousStates) = dlsym(so_hndl,"fmi2GetContinuousStates");
	
	fmi2CallbackFunctions tmp = {fmilogger,calloc,free,NULL,NULL};
	callbackFuncs = new fmi2CallbackFunctions(tmp);
	// Initiate the FMI component using GUID from XML
	m = _fmi2Instantiate(modelname,fmi2ModelExchange,guid,"",callbackFuncs,fmi2False,fmi2False);
}
	
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
    // The number of FMUs
	std::istringstream isss( argv[3] );
        int numOfFMUs;
        if (isss >> numOfFMUs)
        {
            // Conversion successful
        } 
    // The number of Cores
	std::istringstream issss( argv[4] );
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
		// Create FMI object
		vector<FMI*> fmus;
		//#pragma omp parallel for
		for(int i = 0;i<numOfFMUs;i++){
	
			fmus.push_back(new FMI("SerialTanks","{8c4e810f-3df3-4a00-8276-176fa3c9f9e0}",1,0,"binaries/linux64/SerialTanks.so"));
			fmus[i]->fmi_valu = new fmi2Real[numOfTanks];
			fmus[i]->fmi_valx = new fmi2Real[numOfTanks];
			fmus[i]->fmi_vala = new fmi2Real[numOfTanks];
			fmus[i]->fmi_valder = new fmi2Real[numOfTanks];	
		}
		// set the start and end time of simulation
		double Tstart = 0;
		double Tend = 1000.0;
		// fixed step size of 10 milli-seconds
		double dt = 0.1; 
		double time  = Tstart;
		
		// Arrays of inputs, outputs, and variables
	  	// determine continuous and discrete states
	 	fmi2ValueReference *fmi_refu = new fmi2ValueReference[numOfTanks];
		fmi2ValueReference *fmi_refa = new fmi2ValueReference[numOfTanks];
		fmi2ValueReference *fmi_refx = new fmi2ValueReference[numOfTanks];			   
		fmi2ValueReference *fmi_refder = new fmi2ValueReference[numOfTanks];
		
		size_t nx = numOfTanks;
		// Assign value references
		for(int i = 0;i<numOfTanks;i++){
			
			fmi_refu[i]=numOfTanks*2+i;
			fmi_refa[i]=numOfTanks*3+i;
			fmi_refder[i]=numOfTanks+i;
			fmi_refx[i]=i;
		}
		#pragma omp parallel for
		for(int i = 0;i<numOfFMUs;i++){
			// FMU component
			fmi2Component m = fmus[i]->m;
    		// Set Time
			(fmus[i])->_fmi2SetTime(m, time);
			// initialize
   			// determine continuous and discrete states
   			(fmus[i])->_fmi2SetupExperiment(m,fmi2False,0.0, Tstart, fmi2True,Tend);
   			(fmus[i])->_fmi2EnterInitializationMode(m);
   			(fmus[i])->_fmi2ExitInitializationMode(m);
			// Get first values
		  	(fmus[i])->_fmi2GetReal(m,fmi_refu,nx,fmus[i]->fmi_valu);
		   	(fmus[i])->_fmi2GetReal(m,fmi_refa,nx,fmus[i]->fmi_vala);
		   	(fmus[i])->_fmi2GetReal(m,fmi_refx,nx,fmus[i]->fmi_valx);
		   	(fmus[i])->_fmi2GetReal(m,fmi_refder,nx,fmus[i]->fmi_valder);
			// enter Continuous-Time Mode
			(fmus[i])->_fmi2EnterContinuousTimeMode(m);
   		}
	   	while (time < Tend)
	   	{	
	   		// advance time
			time = time + dt;
			// Map inputs and outputs of FMUs
			#pragma omp parallel for
	   		for(int i = 0;i<numOfFMUs;i++){
				
				// FMU component
				fmi2Component m = fmus[i]->m;
				// Not the First FMU - Exception! To be able to couple FMUS
				if (i > 0){
					
					fmi2Real val;	
					// Get input values from the previous FMU
					val = fmus[i-1]->fmi_valx[numOfTanks-1]*fmus[i-1]->fmi_vala[numOfTanks-1];
					fmus[i]->fmi_valu[0] = val;		
					(fmus[i])->_fmi2SetReal(m,fmi_refu,nx,fmus[i]->fmi_valu);
				}
			}
			#pragma omp parallel for
			for(int i = 0;i<numOfFMUs;i++){
				// FMU component
				fmi2Component m = fmus[i]->m;
				// forward Euler method
				// compute derivatives
				(fmus[i])->_fmi2GetDerivatives(m, fmus[i]->fmi_valder,nx);
				for (int j = 0;j<numOfTanks;j++){
					// Exception! To be able to couple FMUS
					if (i > 0 && j == 0){
						
						fmus[i]->fmi_valx[j] = fmus[i]->fmi_valx[j] + dt*(fmus[i]->fmi_valu[j]-fmus[i]->fmi_valx[j]*fmus[i]->fmi_vala[j]);
					}
					else {
						// compute derivatives
						fmus[i]->fmi_valx[j] = fmus[i]->fmi_valx[j] + dt*fmus[i]->fmi_valder[j];
					}
				}		
				(fmus[i])->_fmi2SetContinuousStates(m, fmus[i]->fmi_valder, nx); 
				(fmus[i])->_fmi2SetReal(m,fmi_refx,nx,fmus[i]->fmi_valx);
				// Set time
				(fmus[i])->_fmi2SetTime(m, time);
				// print results
				//cout<<" "<<time<<" FMU"<<i<<" "<<fmus[i]->fmi_valx[0]<<" "<<fmus[i]->fmi_valx[1];
			}
			//cout<<endl;	
		}
		
		// cleanup
		for(int i = 0;i<numOfFMUs;i++){
	
			delete fmus[i];
		}

		double end = omp_get_wtime();
		double seconds = end-start;	 
		cout <<seconds<< endl;
	}	
}

