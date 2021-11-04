#ifndef ParallelComponentsComponentsZone_h_
#define ParallelComponentsComponentsZone_h_
#include "fmi.h"

class ParallelComponentsComponentsZone:
	public FMI
{
	public:
		ParallelComponentsComponentsZone():
			FMI
			(
				"ParallelComponentsComponentsZone",
				"{f68b854b-7124-4b2c-b3fe-e075808655fe}",
				2,
				14,
				"/home/ozi/Documents/parallel_computing/parallelFMI/threeComponentRoom/binaries/linux64/ParallelComponents_Components_Zone.so"
			)
		{
		}
		double get_V() { return get_real(100663296); }
		void set_V(double val) { set_real(100663296,val); }
		double get_QRooInt_flow() { return get_real(16777216); }
		void set_QRooInt_flow(double val) { set_real(16777216,val); }
		double get_QRooC_flow_nominal() { return get_real(100663297); }
		void set_QRooC_flow_nominal(double val) { set_real(100663297,val); }
		double get_mA_flow_nominal() { return get_real(100663298); }
		void set_mA_flow_nominal(double val) { set_real(100663298,val); }
		double get_dTFan() { return get_real(16777217); }
		void set_dTFan(double val) { set_real(16777217,val); }
		double get_QCoiC_flow_nominal() { return get_real(100663299); }
		void set_QCoiC_flow_nominal(double val) { set_real(100663299,val); }
		double get_TASup_nominal() { return get_real(16777218); }
		void set_TASup_nominal(double val) { set_real(16777218,val); }
		double get_TRooSet() { return get_real(16777219); }
		void set_TRooSet(double val) { set_real(16777219,val); }
		double get_TOut_nominal() { return get_real(16777220); }
		void set_TOut_nominal(double val) { set_real(16777220,val); }
		double get_THeaRecLvg() { return get_real(100663300); }
		void set_THeaRecLvg(double val) { set_real(100663300,val); }
		double get_eps() { return get_real(16777221); }
		void set_eps(double val) { set_real(16777221,val); }
		int get_vol_energyDynamics() { return get_int(100663301); }
		void set_vol_energyDynamics(int val) { set_int(100663301,val); }
		int get_vol_massDynamics() { return get_int(100663302); }
		void set_vol_massDynamics(int val) { set_int(100663302,val); }
		int get_vol_substanceDynamics() { return get_int(100663303); }
		void set_vol_substanceDynamics(int val) { set_int(100663303,val); }
		int get_vol_traceDynamics() { return get_int(100663304); }
		void set_vol_traceDynamics(int val) { set_int(100663304,val); }
		double get_vol_p_start() { return get_real(100663305); }
		void set_vol_p_start(double val) { set_real(100663305,val); }
		double get_vol_T_start() { return get_real(16777222); }
		void set_vol_T_start(double val) { set_real(16777222,val); }
		double get_vol_X_start_1_() { return get_real(16777223); }
		void set_vol_X_start_1_(double val) { set_real(16777223,val); }
		double get_vol_X_start_2_() { return get_real(16777224); }
		void set_vol_X_start_2_(double val) { set_real(16777224,val); }
		double get_vol_mSenFac() { return get_real(100663306); }
		void set_vol_mSenFac(double val) { set_real(100663306,val); }
		bool get_vol_prescribedHeatFlowRate() { return get_bool(100663308); }
		void set_vol_prescribedHeatFlowRate(bool val) { set_bool(100663308,val); }
		bool get_vol_simplify_mWat_flow() { return get_bool(100663309); }
		void set_vol_simplify_mWat_flow(bool val) { set_bool(100663309,val); }
		bool get_vol_use_C_flow() { return get_bool(100663310); }
		void set_vol_use_C_flow(bool val) { set_bool(100663310,val); }
		double get_vol_m_flow_nominal() { return get_real(100663311); }
		void set_vol_m_flow_nominal(double val) { set_real(100663311,val); }
		int get_vol_nPorts() { return get_int(100663312); }
		void set_vol_nPorts(int val) { set_int(100663312,val); }
		double get_vol_m_flow_small() { return get_real(100663313); }
		void set_vol_m_flow_small(double val) { set_real(100663313,val); }
		bool get_vol_allowFlowReversal() { return get_bool(100663314); }
		void set_vol_allowFlowReversal(bool val) { set_bool(100663314,val); }
		double get_vol_V() { return get_real(100663315); }
		void set_vol_V(double val) { set_real(100663315,val); }
		double get_vol_ports_1__m_flow() { return get_real(352321536); }
		void set_vol_ports_1__m_flow(double val) { set_real(352321536,val); }
		double get_vol_ports_1__p() { return get_real(100663316); }
		void set_vol_ports_1__p(double val) { set_real(100663316,val); }
		double get_vol_ports_1__h_outflow() { return get_real(369098773); }
		void set_vol_ports_1__h_outflow(double val) { set_real(369098773,val); }
		double get_vol_ports_1__Xi_outflow_1_() { return get_real(335544323); }
		void set_vol_ports_1__Xi_outflow_1_(double val) { set_real(335544323,val); }
		double get_vol_ports_2__p() { return get_real(100663318); }
		void set_vol_ports_2__p(double val) { set_real(100663318,val); }
		double get_vol_ports_2__h_outflow() { return get_real(369098773); }
		void set_vol_ports_2__h_outflow(double val) { set_real(369098773,val); }
		double get_vol_ports_2__Xi_outflow_1_() { return get_real(335544323); }
		void set_vol_ports_2__Xi_outflow_1_(double val) { set_real(335544323,val); }
		double get_vol_heatPort_T() { return get_real(335544320); }
		void set_vol_heatPort_T(double val) { set_real(335544320,val); }
		double get_vol_heatPort_Q_flow() { return get_real(637534231); }
		void set_vol_heatPort_Q_flow(double val) { set_real(637534231,val); }
		double get_vol_T() { return get_real(335544320); }
		void set_vol_T(double val) { set_real(335544320,val); }
		double get_vol_U() { return get_real(33554432); }
		void set_vol_U(double val) { set_real(33554432,val); }
		double get_vol_p() { return get_real(100663320); }
		void set_vol_p(double val) { set_real(100663320,val); }
		double get_vol_m() { return get_real(100663321); }
		void set_vol_m(double val) { set_real(100663321,val); }
		double get_vol_Xi_1_() { return get_real(335544323); }
		void set_vol_Xi_1_(double val) { set_real(335544323,val); }
		double get_vol_mXi_1_() { return get_real(33554433); }
		void set_vol_mXi_1_(double val) { set_real(33554433,val); }
		double get_vol_dynBal_U() { return get_real(33554432); }
		void set_vol_dynBal_U(double val) { set_real(33554432,val); }
		double get_der_vol_dynBal_U_() { return get_real(587202560); }
		void set_der_vol_dynBal_U_(double val) { set_real(587202560,val); }
		double get_vol_dynBal_mXi_1_() { return get_real(33554433); }
		void set_vol_dynBal_mXi_1_(double val) { set_real(33554433,val); }
		double get_der_vol_dynBal_mXi_1__() { return get_real(587202561); }
		void set_der_vol_dynBal_mXi_1__(double val) { set_real(587202561,val); }
		double get_theCon_Q_flow() { return get_real(637534299); }
		void set_theCon_Q_flow(double val) { set_real(637534299,val); }
		double get_theCon_dT() { return get_real(637534300); }
		void set_theCon_dT(double val) { set_real(637534300,val); }
		double get_theCon_port_a_T() { return get_real(16777228); }
		void set_theCon_port_a_T(double val) { set_real(16777228,val); }
		double get_theCon_port_a_Q_flow() { return get_real(637534299); }
		void set_theCon_port_a_Q_flow(double val) { set_real(637534299,val); }
		double get_theCon_port_b_T() { return get_real(335544320); }
		void set_theCon_port_b_T(double val) { set_real(335544320,val); }
		double get_theCon_port_b_Q_flow() { return get_real(637534329); }
		void set_theCon_port_b_Q_flow(double val) { set_real(637534329,val); }
		double get_theCon_G() { return get_real(16777225); }
		void set_theCon_G(double val) { set_real(16777225,val); }
		double get_preHea_Q_flow() { return get_real(100663389); }
		void set_preHea_Q_flow(double val) { set_real(100663389,val); }
		double get_preHea_T_ref() { return get_real(16777226); }
		void set_preHea_T_ref(double val) { set_real(16777226,val); }
		double get_preHea_alpha() { return get_real(16777227); }
		void set_preHea_alpha(double val) { set_real(16777227,val); }
		double get_preHea_port_T() { return get_real(335544320); }
		void set_preHea_port_T(double val) { set_real(335544320,val); }
		double get_preHea_port_Q_flow() { return get_real(637534302); }
		void set_preHea_port_Q_flow(double val) { set_real(637534302,val); }
		double get_senTemRoo_T() { return get_real(335544320); }
		void set_senTemRoo_T(double val) { set_real(335544320,val); }
		double get_senTemRoo_port_T() { return get_real(335544320); }
		void set_senTemRoo_port_T(double val) { set_real(335544320,val); }
		double get_senTemRoo_port_Q_flow() { return get_real(100663391); }
		void set_senTemRoo_port_Q_flow(double val) { set_real(100663391,val); }
		double get_y() { return get_real(335544320); }
		void set_y(double val) { set_real(335544320,val); }
		bool get_outlet1_use_p_in() { return get_bool(100663403); }
		void set_outlet1_use_p_in(bool val) { set_bool(100663403,val); }
		bool get_outlet1_allowFlowReversal() { return get_bool(100663404); }
		void set_outlet1_allowFlowReversal(bool val) { set_bool(100663404,val); }
		double get_outlet1_m_flow() { return get_real(335544321); }
		void set_outlet1_m_flow(double val) { set_real(335544321,val); }
		double get_outlet1_forward_T() { return get_real(335544322); }
		void set_outlet1_forward_T(double val) { set_real(335544322,val); }
		double get_outlet1_forward_X_w() { return get_real(335544323); }
		void set_outlet1_forward_X_w(double val) { set_real(335544323,val); }
		bool get_inlet1_use_p_in() { return get_bool(100663415); }
		void set_inlet1_use_p_in(bool val) { set_bool(100663415,val); }
		bool get_inlet1_allowFlowReversal() { return get_bool(100663416); }
		void set_inlet1_allowFlowReversal(bool val) { set_bool(100663416,val); }
		double get_inlet1_m_flow() { return get_real(352321536); }
		void set_inlet1_m_flow(double val) { set_real(352321536,val); }
		double get_inlet1_forward_T() { return get_real(352321537); }
		void set_inlet1_forward_T(double val) { set_real(352321537,val); }
		double get_inlet1_forward_X_w() { return get_real(352321538); }
		void set_inlet1_forward_X_w(double val) { set_real(352321538,val); }
		double get_fixedTemperature_T() { return get_real(16777228); }
		void set_fixedTemperature_T(double val) { set_real(16777228,val); }
		double get_fixedTemperature_port_T() { return get_real(16777228); }
		void set_fixedTemperature_port_T(double val) { set_real(16777228,val); }
		double get_fixedTemperature_port_Q_flow() { return get_real(637534329); }
		void set_fixedTemperature_port_Q_flow(double val) { set_real(637534329,val); }
};

#endif