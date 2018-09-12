//# define ARMA_DONT_USE_WRAPPER
# include <armadillo>
# include <iostream>
# include "fun_return.h"
#include "load_system_data.h"

using namespace std;
using namespace arma;

VCtrig VoltReg_trigger(mat Dl, cx_mat Z, int Nnum, int timer_overvolt, int timer_undervolt)
{
	VCtrig vc_trig;
	
	sysdata sysinfo = load_system_data();
	int Ldl = sysinfo.Dl.n_rows;

	//mat Dl = sysinfo.Dl;
	//cx_mat Z = sysinfo.Z;

	VPQ dpf_temp = DPF_return123(Dl, Z);
	mat Vpolar = dpf_temp.Vpolar;
	int Lvp = Vpolar.n_rows;

	int j, ja, jb, jc, ia, ib, ic;

	j = 1;
	ja = 0;
	jb = 0;
	jc = 0;

	int cnt_nodes = 0;
	int Lla = 0;
	int Llb = 0;
	int Llc = 0;
	for (int i = 0; i < Ldl; ++i)
	{
		if ((int)Dl(i, 0) != 0)
			cnt_nodes = cnt_nodes + 1;
		if ((int)Dl(i, 6) != 0)
			Lla = Lla + 1;
		if ((int)Dl(i, 8) != 0)
			Llb = Llb + 1;
		if ((int)Dl(i, 10) != 0)
			Llc = Llc + 1;
	}

	cnt_nodes = cnt_nodes + 1;//No of nodes= No of branches +1
	mat Node_f = zeros(1, cnt_nodes);
	mat Load_a = zeros(1, Lla);
	mat Load_b = zeros(1, Llb);
	mat Load_c = zeros(1, Llc);

	for (int i = 0; i < Ldl && j<cnt_nodes && ja<Lla && jb<Llb &&jc<Llc; ++i)
	{
		if ((int)Dl(i, 2) != 0)
		{
			Node_f(0, j) = Dl(i, 2);
			++j;
		}
		if ((int)Dl(i, 6) != 0)
		{
			Load_a(ja) = Dl(i, 2);
			++ja;
		}
		if ((int)Dl(i, 8) != 0)
		{
			Load_b(jb) = Dl(i, 2);
			++jb;
		}
		if ((int)Dl(i, 10) != 0)
		{
			Load_c(jc) = Dl(i, 2);
			++jc;
		}
	}
	Node_f = Node_f.st();

	y_re Y_return = form_Y_abc(Dl, sysinfo.Z, sysinfo.bkva, sysinfo.bkv);

	int Lnum_a = Y_return.Lnum_a;
	int Lnum_b = Y_return.Lnum_b;
	int Lnum_c = Y_return.Lnum_c;

	Vabc V_abc = V_abc_list(Vpolar, Node_f, Lvp, Lnum_a, Lnum_b, Lnum_c);

	int Lna = V_abc.Lna;
	int Lnb = V_abc.Lnb;
	int Lnc = V_abc.Lnc;
	mat Node_a = V_abc.Node_a;
	mat Node_b = V_abc.Node_b;
	mat Node_c = V_abc.Node_c;
	mat V_a = V_abc.V_a;
	mat V_b = V_abc.V_b;
	mat V_c = V_abc.V_c;
	mat theta_a = V_abc.theta_a;
	mat theta_b = V_abc.theta_b;
	mat theta_c = V_abc.theta_c;

	mat Vmin_abc, Vmax_abc;
	Vmin_abc << min(min(V_a)) << min(min(V_b)) << min(min(V_c)) << endr;
	Vmax_abc << max(max(V_a)) << max(max(V_b)) << max(max(V_c)) << endr;
	double Vmin_orig = min(min(Vmin_abc));
	double Vmax_orig = max(max(Vmax_abc));
	//cout << "Vmax (p.u.) = " << Vmax_orig << endl;
	//cout << "Vmin (p.u.) = " << Vmin_orig << endl;

	int timer_overv, timer_underv, re, trigger;
	if (Vmax_orig > 1.05 + 0.0001)
	{
		re = 1;
		timer_overv = timer_overvolt + 1;
		timer_underv = 0;
	}
	else if (Vmin_orig < 0.95 - 0.0001)
	{
		re = -1;
		timer_underv = timer_undervolt + 1;
		timer_overv = 0;
	}
	else
	{
		re = 0;
		timer_overv = 0;
		timer_underv = 0;
	}

	if (timer_overv == 3 || timer_underv == 3)
	//if (timer_overv == 1 || timer_underv == 1) //For checking
	{
		trigger = 1;
		timer_overv = 0;
		timer_underv = 0;
	}
	else
	{
		trigger = 0;
	}

	vc_trig.re1 = re;
	vc_trig.trigger = trigger;
	vc_trig.timer_overv = timer_overv;
	vc_trig.timer_underv = timer_underv;
	return vc_trig;
}