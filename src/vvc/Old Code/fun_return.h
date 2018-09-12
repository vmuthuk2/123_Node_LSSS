#ifndef FUN_RETURN_HPP_
#define FUN_RETURN_HPP_
# define ARMA_DONT_USE_WRAPPER

#include <armadillo>
//using namespace arma;
struct y_re   //output of form_Yabc_34(*)
{
	arma::cx_mat Y_a;
	arma::cx_mat Y_b;
	arma::cx_mat Y_c;
	arma::cx_mat brnches;
	int Nnum;
	arma::mat Lnum;
	int Lnum_a;
	int Lnum_b;
	int Lnum_c;
};

//y_re form_Y_abc_34(mat Dl, cx_mat Z, double bkva, double bkv);//zhan's
y_re form_Y_abc(arma::mat Dl, arma::cx_mat Z, double bkva, double bkv);//Yue's

struct Vabc // document the valid Node Voltages and corrsponding node number for each phase
{
	arma::mat V_a, V_b, V_c;
	arma::mat theta_a, theta_b, theta_c;
	arma::mat Node_a, Node_b, Node_c;
	int Lna, Lnb, Lnc;
};

Vabc V_abc_list(arma::mat Vpolar, arma::mat Node_f,int Lvp, int Lnum_a, int Lnum_b, int Lnum_c);

struct newbrn	//output of rename_brn(*)
{
	arma::cx_mat newbrn_a;
	arma::cx_mat newbrn_b;
	arma::cx_mat newbrn_c;

};

newbrn rename_brn(arma::mat Node_a, arma::mat Node_b, arma::mat Node_c, arma::cx_mat brn_a, arma::cx_mat brn_b, arma::cx_mat brn_c, int Lnum_a, int Lnum_b, int Lnum_c,int Lna, int Lnb, int Lnc );


struct VPQ //output of DPF ( Voltage PQb and PQL)
{
	arma::mat Vpolar;
	arma::mat PQb;
	arma::mat PQL;
	arma::mat Ib;
	arma::mat IL;
	arma::mat Qset_a,Qset_b,Qset_c;
};

VPQ DPF_return123(arma::mat Dl, arma::cx_mat Z);//added by Valli


arma::mat form_Ftheta(arma::cx_mat Y, arma::mat V, arma::mat theta, arma::cx_mat brn, int Ln, int Lnm);//output of form_Ftheta(*)
arma::mat form_Fv(arma::cx_mat Y, arma::mat V, arma::mat theta, arma::cx_mat brn, int Ln, int Lnm);//output of form_Fv(*)
arma::mat form_J(arma::cx_mat Y, arma::mat V, arma::mat theta, int Lnm);//out put of form_J(*)
arma::mat form_L(arma::cx_mat Y, arma::mat V, arma::mat theta, int Lnm);//out put of form_L(*)

struct VCtrig //output of VoltReg_trigger function
{
	int trigger;
	int timer_overv;
	int timer_underv;
	int re1;
};

VCtrig VoltReg_trigger(arma::mat Dl, arma::cx_mat Z, int Nnum, int timer_overvolt, int timer_undervolt);//output of VoltReg_trigger

struct VRreg //output of VoltageRegulation_v4
{
	arma::mat Dl_vr;
	arma::mat re_vr;
};

VRreg volt_reg_solver(arma::mat Dl, int Nnum, double v0);//output of volt_reg_solver

struct vrmove
{
	arma::mat VR1;
	arma::mat VR2;
	arma::mat VR3;
	arma::mat VR4;
};

vrmove VR_move(arma::mat Vpolar, arma::mat Dl, arma::mat v_lim, int Lvp);

struct tapsearch
{
	arma::mat aVR1;
	arma::mat aVR2;
	arma::mat aVR3;
	arma::mat aVR4;
	arma::mat results;
};

tapsearch tap_searching_first(arma::mat VR1, arma::mat VR2, arma::mat VR3, arma::mat VR4, arma::mat Dl_bf_voltc, double v0, int Nnum, arma::cx_mat Z, arma::mat Node_f, int Lnum_a, int Lnum_b, int Lnum_c);

struct gradV
{
	double Vmin;
	double Vmax;
	arma::mat Dl_volt;
	double Ploss_afterVR;
	int trigger1;
	int no_itea;
};

gradV grad_Volt_opt_solve(arma::mat Dl_vr, double v0, arma::cx_mat Z, double bkva);//output of grad_Volt_opt_solve

//struct qadj
//{
//	arma::mat Dl_adj;
//	arma::mat flag;
//};
//
//qadj qinj_adjustment(arma::mat Dl_check);

struct gradVAR
{
	double V_mn;
	double V_mx;
	double totalloss;
	arma::mat Qcmd;
	int itea_stepsch;
};

gradVAR var_opt_loss_solver(arma::mat Dl_vr, double v0, arma::cx_mat Z, double bkva);//output of var_opt_loss_solver
#endif