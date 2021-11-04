#ifndef ParallelComponentsComponentsController_h_
#define ParallelComponentsComponentsController_h_
#include "fmi.h"

class ParallelComponentsComponentsController:
	public FMI
{
	public:
		ParallelComponentsComponentsController():
			FMI
			(
				"ParallelComponentsComponentsController",
				"{c41e4bef-4e58-404e-845f-80f75bd76639}",
				0,
				4,
				"/home/ozi/Documents/parallel_computing/parallelFMI/threeComponentRoom/binaries/linux64/ParallelComponents_Components_Controller.so"
			)
		{
		}
		double get_eps() { return get_real(16777216); }
		void set_eps(double val) { set_real(16777216,val); }
		double get_TASup_nominal() { return get_real(16777217); }
		void set_TASup_nominal(double val) { set_real(16777217,val); }
		double get_TRooSet() { return get_real(16777218); }
		void set_TRooSet(double val) { set_real(16777218,val); }
		double get_TOut_nominal() { return get_real(16777219); }
		void set_TOut_nominal(double val) { set_real(16777219,val); }
		double get_THeaRecLvg() { return get_real(100663296); }
		void set_THeaRecLvg(double val) { set_real(100663296,val); }
		double get_QRooInt_flow() { return get_real(16777220); }
		void set_QRooInt_flow(double val) { set_real(16777220,val); }
		double get_QRooC_flow_nominal() { return get_real(100663297); }
		void set_QRooC_flow_nominal(double val) { set_real(100663297,val); }
		double get_mA_flow_nominal() { return get_real(100663298); }
		void set_mA_flow_nominal(double val) { set_real(100663298,val); }
		double get_dTFan() { return get_real(16777221); }
		void set_dTFan(double val) { set_real(16777221,val); }
		double get_QCoiC_flow_nominal() { return get_real(100663299); }
		void set_QCoiC_flow_nominal(double val) { set_real(100663299,val); }
		double get_TWSup_nominal() { return get_real(16777222); }
		void set_TWSup_nominal(double val) { set_real(16777222,val); }
		double get_TWRet_nominal() { return get_real(16777223); }
		void set_TWRet_nominal(double val) { set_real(16777223,val); }
		double get_mW_flow_nominal() { return get_real(100663300); }
		void set_mW_flow_nominal(double val) { set_real(100663300,val); }
		double get_con_reference() { return get_real(100663302); }
		void set_con_reference(double val) { set_real(100663302,val); }
		double get_con_u() { return get_real(352321536); }
		void set_con_u(double val) { set_real(352321536,val); }
		bool get_con_y() { return get_bool(369098757); }
		void set_con_y(bool val) { set_bool(369098757,val); }
		double get_con_bandwidth() { return get_real(16777224); }
		void set_con_bandwidth(double val) { set_real(16777224,val); }
		bool get_con_pre_y_start() { return get_bool(16777225); }
		void set_con_pre_y_start(bool val) { set_bool(16777225,val); }
		double get_TRooSetPoi_k() { return get_real(100663302); }
		void set_TRooSetPoi_k(double val) { set_real(100663302,val); }
		double get_TRooSetPoi_y() { return get_real(100663302); }
		void set_TRooSetPoi_y(double val) { set_real(100663302,val); }
		bool get_mWat_flow_u() { return get_bool(369098757); }
		void set_mWat_flow_u(bool val) { set_bool(369098757,val); }
		double get_mWat_flow_realTrue() { return get_real(16777226); }
		void set_mWat_flow_realTrue(double val) { set_real(16777226,val); }
		double get_mWat_flow_realFalse() { return get_real(100663303); }
		void set_mWat_flow_realFalse(double val) { set_real(100663303,val); }
		double get_mWat_flow_y() { return get_real(335544320); }
		void set_mWat_flow_y(double val) { set_real(335544320,val); }
		double get_u() { return get_real(352321536); }
		void set_u(double val) { set_real(352321536,val); }
		double get_y() { return get_real(335544320); }
		void set_y(double val) { set_real(335544320,val); }
};

#endif