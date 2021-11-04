#ifndef _fmi_h_
#define _fmi_h_
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include "time.h"
#include <cstdlib>
#include <stdlib.h>
#include <dlfcn.h>
#include "assert.h"
#include <vector>
#include "fmi2Functions.h"
#include "fmi2FunctionTypes.h"
#include "fmi2TypesPlatform.h"

using namespace std;

class FMI {

public:
		/**
		 * This constructs a wrapper around an FMI. The constructor
		 * must be provided with the FMI's GUID, the number of state variables,
		 * number of event indicators, and the path to the .so file
		 * that contains the FMI functions for this model.
		 */	
		// Get the value of a real variable
		double get_real(int k);
		// Set the value of a real variable
		void set_real(int k, double val);
		// Get the value of an integer variable
		int get_int(int k);
		// Set the value of an integer variable
		void set_int(int k, int val);
		// Get the value of a boolean variable
		bool get_bool(int k);
		// Set the value of a boolean variable
		void set_bool(int k, bool val);
		/// Event iterations
		virtual void do_events(double* q, double Tend, double Tnext);	
		/// Copy the initial state of the model to q
		virtual void init(double* q, double Tstart, double Tend);
		/// Compute the derivative for state q and put it in dq
		virtual void der_func(const double* q, double* dq, double time);
		/// Destructor
		virtual ~FMI();
		/// Iterate events
		void iterate_events();
		// Create the FMI component
		fmi2CallbackFunctions tmp = {FMI::fmilogger,calloc,free,NULL,NULL};
		fmi2CallbackFunctions* callbackFuncs;
		fmi2Component c;
			
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
		fmi2Status (*_fmi2Terminate)(fmi2Component);
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

		// Instant of the next time event
		double next_time_event;
		// so library handle
		void* so_hndl;
		// Are we in continuous time mode?
		bool cont_time_mode;
		// Number of event indicators that are not governed by the FMI
		int num_extra_event_indicators;

		fmi2Real* q;
		fmi2Real* dq;
		fmi2Real* z;
		int numStates;
		int numEvents;

		static void fmilogger(
			fmi2ComponentEnvironment componentEnvironment,
			fmi2String instanceName,
			fmi2Status status,
			fmi2String category,
			fmi2String message, ...)
		{
			std::cerr << message << std::endl;
		}	

};

FMI::FMI(const char* modelname,const char* guid,int num_state_variables,int num_event_indicators,const char* so_file_name,const double tolerance)
{
	so_hndl = dlopen(so_file_name, RTLD_LAZY);
	if (!so_hndl){
		cout<<"Could not load so file"<< endl;
    	}
	
	//*(void**) (&_fmi2Instantiate) = dlsym(so_hndl,"fmi2Instantiate");
	// This only works with a POSIX compliant compiler/system
	_fmi2Instantiate = (fmi2Component (*)(fmi2String, fmi2Type,
		fmi2String, fmi2String, const fmi2CallbackFunctions*,
		fmi2Boolean, fmi2Boolean))dlsym(so_hndl,"fmi2Instantiate");
	assert(_fmi2Instantiate != NULL);
	_fmi2FreeInstance = (void (*)(fmi2Component))dlsym(so_hndl,"fmi2FreeInstance");
	assert(_fmi2FreeInstance != NULL);
	_fmi2SetupExperiment = (fmi2Status (*)(fmi2Component, fmi2Boolean,
		fmi2Real, fmi2Real, fmi2Boolean, fmi2Real))dlsym(so_hndl,"fmi2SetupExperiment");
	assert(_fmi2SetupExperiment != NULL);
	_fmi2EnterInitializationMode = (fmi2Status (*)(fmi2Component))dlsym(so_hndl,"fmi2EnterInitializationMode");
	assert(_fmi2EnterInitializationMode != NULL);
	_fmi2ExitInitializationMode = (fmi2Status (*)(fmi2Component))dlsym(so_hndl,"fmi2ExitInitializationMode");
	assert(_fmi2ExitInitializationMode != NULL);
	_fmi2GetReal = (fmi2Status (*)(fmi2Component, const fmi2ValueReference*, size_t, fmi2Real*))
		dlsym(so_hndl,"fmi2GetReal");
	assert(_fmi2GetReal != NULL);
	_fmi2GetInteger = (fmi2Status (*)(fmi2Component, const fmi2ValueReference*, size_t, fmi2Integer*)) 
		dlsym(so_hndl,"fmi2GetInteger");
	assert(_fmi2GetInteger != NULL);
	_fmi2GetBoolean = (fmi2Status (*)(fmi2Component, const fmi2ValueReference*, size_t, fmi2Boolean*))
		dlsym(so_hndl,"fmi2GetBoolean");
	assert(_fmi2GetBoolean != NULL);
	_fmi2GetString = (fmi2Status (*)(fmi2Component, const fmi2ValueReference*, size_t, fmi2String*))
		dlsym(so_hndl,"fmi2GetString");
	assert(_fmi2GetString != NULL);
	_fmi2SetReal = (fmi2Status (*)(fmi2Component, const fmi2ValueReference*, size_t, const fmi2Real*))
		dlsym(so_hndl,"fmi2SetReal");
	assert(_fmi2SetReal != NULL);
	_fmi2SetInteger = (fmi2Status (*)(fmi2Component, const fmi2ValueReference*, size_t, const fmi2Integer*))
		dlsym(so_hndl,"fmi2SetInteger");
	assert(_fmi2SetInteger != NULL);
	_fmi2SetBoolean = (fmi2Status (*)(fmi2Component, const fmi2ValueReference*, size_t, const fmi2Boolean*))
		dlsym(so_hndl,"fmi2SetBoolean");
	assert(_fmi2SetBoolean != NULL);
	_fmi2SetString = (fmi2Status (*)(fmi2Component, const fmi2ValueReference*, size_t, const fmi2String*))
		dlsym(so_hndl,"fmi2SetString");
	assert(_fmi2SetString != NULL);
	_fmi2EnterEventMode = (fmi2Status (*)(fmi2Component))dlsym(so_hndl,"fmi2EnterEventMode");
	assert(_fmi2EnterEventMode != NULL);
	_fmi2NewDiscreteStates = (fmi2Status (*)(fmi2Component,fmi2EventInfo*))dlsym(so_hndl,"fmi2NewDiscreteStates");
	assert(_fmi2NewDiscreteStates != NULL);
	_fmi2EnterContinuousTimeMode = (fmi2Status (*)(fmi2Component))dlsym(so_hndl,"fmi2EnterContinuousTimeMode");
	assert(_fmi2EnterContinuousTimeMode != NULL);
	_fmi2CompletedIntegratorStep = (fmi2Status (*)(fmi2Component, fmi2Boolean, fmi2Boolean*, fmi2Boolean*))
		dlsym(so_hndl,"fmi2CompletedIntegratorStep");
	assert(_fmi2CompletedIntegratorStep != NULL);
	_fmi2SetTime = (fmi2Status (*)(fmi2Component, fmi2Real))dlsym(so_hndl,"fmi2SetTime");
	assert(_fmi2SetTime != NULL);
	_fmi2SetContinuousStates = (fmi2Status (*)(fmi2Component, const fmi2Real*, size_t))
		dlsym(so_hndl,"fmi2SetContinuousStates");
	assert(_fmi2SetContinuousStates != NULL);
	_fmi2GetDerivatives = (fmi2Status (*)(fmi2Component, fmi2Real*, size_t))dlsym(so_hndl,"fmi2GetDerivatives");
	assert(_fmi2GetDerivatives != NULL);
	_fmi2GetEventIndicators = (fmi2Status (*)(fmi2Component, fmi2Real*, size_t))dlsym(so_hndl,"fmi2GetEventIndicators");
	assert(_fmi2GetEventIndicators != NULL);
	_fmi2GetContinuousStates = (fmi2Status (*)(fmi2Component, fmi2Real*, size_t))dlsym(so_hndl,"fmi2GetContinuousStates");
	assert(_fmi2GetContinuousStates != NULL);
	_fmi2Terminate = (fmi2Status (*)(fmi2Component))dlsym(so_hndl,"fmi2Terminate");
	assert(_fmi2Terminate != NULL);
	next_time_event = 1E+8;
	//so_hndl = NULL;
	cont_time_mode =false;
	fmi2CallbackFunctions tmp = {fmilogger,calloc,free,NULL,NULL};
	callbackFuncs = new fmi2CallbackFunctions(tmp);
	// Initiate the FMI component
	c = _fmi2Instantiate(modelname,fmi2ModelExchange,guid,"",callbackFuncs,fmi2False,fmi2False);
	numStates = num_state_variables;
	numEvents = num_event_indicators;
}

void FMI::iterate_events()
{
	fmi2Status status;
	if (cont_time_mode)
	{
		status = _fmi2EnterEventMode(c);
		assert(status == fmi2OK);
		cont_time_mode = false;
	}
	// Put into consistent initial state
	fmi2EventInfo eventInfo;
	do
	{
		status = _fmi2NewDiscreteStates(c,&eventInfo);
		assert(status == fmi2OK);
	}
	while (eventInfo.newDiscreteStatesNeeded == fmi2True);
	if (eventInfo.nextEventTimeDefined == fmi2True)
		next_time_event = eventInfo.nextEventTime;
	assert(status == fmi2OK);
}

void FMI::init(double* q,double Tstart,double Tend)
{
	fmi2Status status;
	status = _fmi2SetupExperiment(c,fmi2False,0.0,Tstart,fmi2True,Tend);
	// Initialize all variables
	status = _fmi2EnterInitializationMode(c);
	assert(status == fmi2OK);
	// Done with initialization
	status = _fmi2ExitInitializationMode(c);
	assert(status == fmi2OK);
	// Put into consistent initial state
	iterate_events();
	// Enter continuous time mode to start integration
	status = _fmi2EnterContinuousTimeMode(c);
	assert(status == fmi2OK);
	status = _fmi2GetContinuousStates(c,q,numStates);
	assert(status == fmi2OK);
	cont_time_mode = true;
}

void FMI::der_func(const double* q, double* dq, double time)
{
	fmi2Status status;
	if (!cont_time_mode)
	{
		status = _fmi2EnterContinuousTimeMode(c);
		assert(status == fmi2OK);
		cont_time_mode = true;
	}
	//status = _fmi2SetContinuousStates(c,q,numStates);
	//assert(status == fmi2OK);
	status = _fmi2GetDerivatives(c,dq,numStates);
	assert(status == fmi2OK);
}

void FMI::do_events(double* q, double Tend, double Tnext)
{
	fmi2Status status;
	if (cont_time_mode)
	{
		status = _fmi2EnterEventMode(c);
		assert(status == fmi2OK);
		cont_time_mode = false;
	}
	// Process events
	//fmi2EventInfo eventInfo;			
	//fmi2Real* ind = new fmi2Real[numStates];
	//ind  = z;
	status = _fmi2GetEventIndicators(c,z,numEvents);
	assert(status == fmi2OK);
	//bool events = false;
	//for (int t = 0; t<numEvents;t++){
				
	//	if ((z[t]<0 && ind[t]>=0)||(z[t]>=0 && ind[t]<0)){
	//		events = true;
	//		break;
	//	}	
	//}
	fmi2Boolean enterEventMode;
	fmi2Boolean terminateSimulation;
	status = _fmi2CompletedIntegratorStep(c,fmi2True,&enterEventMode,&terminateSimulation);
	assert(status == fmi2OK);
	//if (enterEventMode == fmi2True || events == true){
		
		//do
		//{
		//	status = _fmi2NewDiscreteStates(c,&eventInfo);
		//	assert(status == fmi2OK);
		//	if (eventInfo.terminateSimulation == fmi2True)
		//		_fmi2Terminate(c);
		//	if (eventInfo.nextEventTimeDefined == fmi2True)
		//	Tnext = std::min(eventInfo.nextEventTime, Tend);
		//	else
		//	Tnext = Tend;
		//}
		//while (eventInfo.newDiscreteStatesNeeded == fmi2True);	
	//}
	iterate_events();
	// Update the state variable array
	status = _fmi2GetContinuousStates(c,q,numStates);
	assert(status == fmi2OK);

}

double FMI::get_real(int k)
{
	const fmi2ValueReference ref = k;
	fmi2Real val;
	fmi2Status status = _fmi2GetReal(c,&ref,1,&val);
	assert(status == fmi2OK);
	return val;
}


void FMI::set_real(int k, double val)
{
	const fmi2ValueReference ref = k;
	fmi2Real fmi_val = val;
	fmi2Status status = _fmi2SetReal(c,&ref,1,&fmi_val);
	assert(status == fmi2OK);
}

int FMI::get_int(int k)
{
	const fmi2ValueReference ref = k;
	fmi2Integer val;
	fmi2Status status = _fmi2GetInteger(c,&ref,1,&val);
	assert(status == fmi2OK);
	return val;
}

void FMI::set_int(int k, int val)
{
	const fmi2ValueReference ref = k;
	fmi2Integer fmi_val = val;
	fmi2Status status = _fmi2SetInteger(c,&ref,1,&fmi_val);
	assert(status == fmi2OK);
}

bool FMI::get_bool(int k)
{
	const fmi2ValueReference ref = k;
	fmi2Boolean val;
	fmi2Status status = _fmi2GetBoolean(c,&ref,1,&val);
	assert(status == fmi2OK);
	return (val == fmi2True);
}

void FMI::set_bool(int k, bool val)
{
	const fmi2ValueReference ref = k;
	fmi2Boolean fmi_val = fmi2False;
	if (val) fmi_val = fmi2True;
	fmi2Status status = _fmi2SetBoolean(c,&ref,1,&fmi_val);
	assert(status == fmi2OK);
}

FMI::~FMI(){
		
	_fmi2FreeInstance(c);
	delete callbackFuncs;
	dlclose(so_hndl);
}

#endif
