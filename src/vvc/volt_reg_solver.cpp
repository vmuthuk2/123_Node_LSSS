//# define ARMA_DONT_USE_WRAPPER
# include <armadillo>
# include <iostream>
# include "fun_return.h"
#include "load_system_data.h"

using namespace std;
using namespace arma;

VRreg volt_reg_solver(mat Dl, int Nnum, double v0)
{
	VRreg vr_reg;

	sysdata sysinfo = load_system_data();
	int Ldl = sysinfo.Dl.n_rows;
	int Wdl = sysinfo.Dl.n_cols;

	//mat Dl = sysinfo.Dl;
	cx_mat Z = sysinfo.Z;

	mat lcd_vr1, lcd_vr2, lcd_vr3, lcd_vr4;
	lcd_vr1 << 1 << 1 << 1 << endr;
	lcd_vr2 << 1 << 1 << 1 << endr;
	lcd_vr3 << 1 << 0 << 0 << endr;
	lcd_vr4 << 1 << 0 << 1 << endr;

	//field<mat> vr_lcd(4,1);
	//vr_lcd << lcd_vr1 << lcd_vr2 << lcd_vr3 << lcd_vr4;

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
		//if ((int)Dl(i, 6) != 0)
		//{
		//	Load_a(ja) = Dl(i, 2);
		//	++ja;
		//}
		//if ((int)Dl(i, 8) != 0)
		//{
		//	Load_b(jb) = Dl(i, 2);
		//	++jb;
		//}
		//if ((int)Dl(i, 10) != 0)
		//{
		//	Load_c(jc) = Dl(i, 2);
		//	++jc;
		//}
	}
	Node_f = Node_f.st();

	y_re Y_return = form_Y_abc(Dl, sysinfo.Z, sysinfo.bkva, sysinfo.bkv);
	int Lbr = Y_return.brnches.n_rows;
	int Wbr = Y_return.brnches.n_cols;

	VPQ dpf_temp = DPF_return123(Dl, Z);
	mat Vpolar = dpf_temp.Vpolar;
	mat PQb = dpf_temp.PQb;
	mat PQL = dpf_temp.PQL;

	int Lvp = Vpolar.n_rows;
	int Wvp = Vpolar.n_cols;
	int Lpq = PQb.n_rows;
	int Wpq = PQb.n_cols;
	double Pla_t = sum(PQL.col(0));
	double Plb_t = sum(PQL.col(2));
	double Plc_t = sum(PQL.col(4));

	// calculate power loss
	mat Pload_total, Ploss_orig_ph;
	double Ploss_orig;
	Pload_total << Pla_t << Plb_t << Plc_t << endr;
	Ploss_orig_ph << (PQb(0, 0) - Pla_t) << (PQb(0, 2) - Plb_t) << (PQb(0, 4) - Plc_t) << endr;
	Ploss_orig = accu(Ploss_orig_ph);

	int Lnum_a = Y_return.Lnum_a;
	int Lnum_b = Y_return.Lnum_b;
	int Lnum_c = Y_return.Lnum_c;

	Vabc V_abc = V_abc_list(Vpolar, Node_f, Lvp, Lnum_a, Lnum_b, Lnum_c);

	mat V_a = V_abc.V_a;
	mat V_b = V_abc.V_b;
	mat V_c = V_abc.V_c;

	mat Vmin_abc, Vmax_abc;
	Vmin_abc << min(min(V_a)) << min(min(V_b)) << min(min(V_c)) << endr;
	Vmax_abc << max(max(V_a)) << max(max(V_b)) << max(max(V_c)) << endr;
	double Vmin_orig = min(min(Vmin_abc));
	double Vmax_orig = max(max(Vmax_abc));

	mat v_lim;
	v_lim << 1.05 << 0.95;

	vrmove vr_move = VR_move(Vpolar, Dl, v_lim, Lvp);
	
	mat VR1 = vr_move.VR1;
	mat VR2 = vr_move.VR2;
	mat VR3 = vr_move.VR3;
	mat VR4 = vr_move.VR4;

	mat Dl_bf_voltc = Dl;
	
	int q, p, n, m;
	q = p = n = m = 1;

	tapsearch tap_search = tap_searching_first(VR1, VR2, VR3, VR4, Dl_bf_voltc, v0, Nnum, Z, Node_f, Lnum_a, Lnum_b, Lnum_c);
	cout << "Came out of tapsearching!" << endl;
	mat aVR1, aVR2, aVR3, aVR4, results;
	aVR1 = tap_search.aVR1;
	aVR2 = tap_search.aVR2;
	aVR3 = tap_search.aVR3;
	aVR4 = tap_search.aVR4;
	//system("pause");
	results = tap_search.results;
	//system("pause");
	mat Dl_vr = Dl;
	Dl_vr(0, span(13, 15)) = aVR1;
	Dl_vr(12, span(13, 15)) = aVR2;
	Dl_vr(34, span(13, 15)) = aVR3;
	Dl_vr(153, span(13, 15)) = aVR4;
	
	double qmax, pmax, nmax, mmax;
	//qmax = min(min(((1.1 - aVR1) % vr_lcd( 0 )) / (0.1 / 16)));
	//pmax = min(min(((1.1 - aVR2) % vr_lcd( 1 )) / (0.1 / 16)));
	//nmax = min(min(((1.1 - aVR3) % vr_lcd( 2 )) / (0.1 / 16)));
	//mmax = min(min(((1.1 - aVR4) % vr_lcd( 3 )) / (0.1 / 16)));
	qmax = min(min(((1.1 - aVR1) % lcd_vr1) / (0.1 / 16)));
	pmax = min(min(((1.1 - aVR2) % lcd_vr2) / (0.1 / 16)));
	nmax = min(min(((1.1 - aVR3) % lcd_vr3) / (0.1 / 16)));
	mmax = min(min(((1.1 - aVR4) % lcd_vr4) / (0.1 / 16)));

	double temp = results(0, 0);
	cout << "2 * pow(10,6) = " << 2 * pow(10, 6) << endl;
	while (temp == 2 * pow(10, 6))
	{
		mat Dl_new = Dl_bf_voltc;
		
		if (q <= qmax)
		{
			Dl_new(0, span(13, 15)) = aVR1 + q*(0.1 / 16)*VR1;
			q = q + 1;
		}
		else
		{
			Dl_new(0, span(13, 15)) = aVR1 + qmax*(0.1 / 16)*VR1;
		}

		if (p <= pmax)
		{
			Dl_new(12, span(13, 15)) = aVR2 + p*(0.1 / 16)*VR2;
			p = p + 1;
		}
		else
		{
			Dl_new(12, span(13, 15)) = aVR2 + pmax*(0.1 / 16)*VR2;
		}

		if (n <= nmax)
		{
			Dl_new(34, span(13,15)) = aVR3 + n*(0.1 / 16)*VR3;
			n = n + 1;
		}
		else
		{
			Dl_new(34, span(13, 15)) = aVR3 + nmax*(0.1 / 16)*VR3;
		}

		if (m <= mmax)
		{
			Dl_new(153, span(13, 15)) = aVR4 + m*(0.1 / 16)*VR4;
			m = m + 1;
		}
		else
		{
			Dl_new(153, span(13, 15)) = aVR4 + mmax*(0.1 / 16)*VR4;
		}

		VPQ dpf_temp1 = DPF_return123(Dl_new, Z);
		mat Vpolar1 = dpf_temp1.Vpolar;

		vrmove vr_move1 = VR_move(Vpolar1, Dl_new, v_lim, Lvp);
		VR1 = vr_move1.VR1;
		VR2 = vr_move1.VR2;
		VR3 = vr_move1.VR3;
		VR4 = vr_move1.VR4;

		tapsearch tap_search1 = tap_searching_first(VR1, VR2, VR3, VR4, Dl_new, v0, Nnum, Z, Node_f, Lnum_a, Lnum_b, Lnum_c);
		
		mat results1;
		aVR1 = tap_search1.aVR1;
		aVR2 = tap_search1.aVR2;
		aVR3 = tap_search1.aVR3;
		aVR4 = tap_search1.aVR4;
		results1 = tap_search1.results;
		temp = results1(0, 0);
		cout << "results: " << temp << endl;
		Dl_vr(0, span(13, 15)) = aVR1;
		Dl_vr(12, span(13, 15)) = aVR2;
		Dl_vr(34, span(13, 15)) = aVR3;
		Dl_vr(153, span(13, 15)) = aVR4;
	}
	cout << "Came out of 2e6 loop \n"; 
	//system("pause");
	q = p = n = m = 1;
	//qmax = min(min(((aVR1 - 0.9) % vr_lcd(0)) / (0.1 / 16)));
	//pmax = min(min(((aVR2 - 0.9) % vr_lcd(1)) / (0.1 / 16)));
	//nmax = min(min(((aVR3 - 0.9) % vr_lcd(2)) / (0.1 / 16)));
	//mmax = min(min(((aVR4 - 0.9) % vr_lcd(3)) / (0.1 / 16)));
	qmax = min(min(((aVR1 - 0.9) % lcd_vr1) / (0.1 / 16)));
	pmax = min(min(((aVR2 - 0.9) % lcd_vr2) / (0.1 / 16)));
	nmax = min(min(((aVR3 - 0.9) % lcd_vr3) / (0.1 / 16)));
	mmax = min(min(((aVR4 - 0.9) % lcd_vr4) / (0.1 / 16)));

	while (temp == 1 * pow(10, 6))
	{
		mat Dl_new = Dl_bf_voltc;

		if (q <= qmax)
		{
			Dl_new(0, span(13, 15)) = aVR1 + q*(0.1 / 16)*VR1;
			q = q + 1;
		}
		else
		{
			Dl_new(0, span(13, 15)) = aVR1 + qmax*(0.1 / 16)*VR1;
		}

		if (p <= pmax)
		{
			Dl_new(12, span(13, 15)) = aVR2 + p*(0.1 / 16)*VR2;
			p = p + 1;
		}
		else
		{
			Dl_new(12, span(13, 15)) = aVR2 + pmax*(0.1 / 16)*VR2;
		}

		if (n <= nmax)
		{
			Dl_new(34, span(13, 15)) = aVR3 + n*(0.1 / 16)*VR3;
			n = n + 1;
		}
		else
		{
			Dl_new(34, span(13, 15)) = aVR3 + nmax*(0.1 / 16)*VR3;
		}

		if (m <= mmax)
		{
			Dl_new(153, span(13, 15)) = aVR4 + m*(0.1 / 16)*VR4;
			m = m + 1;
		}
		else
		{
			Dl_new(153, span(13, 15)) = aVR4 + mmax*(0.1 / 16)*VR4;
		}

		VPQ dpf_temp1 = DPF_return123(Dl_new, Z);
		mat Vpolar1 = dpf_temp1.Vpolar;

		vrmove vr_move1 = VR_move(Vpolar1, Dl_new, v_lim, Lvp);
		VR1 = vr_move1.VR1;
		VR2 = vr_move1.VR2;
		VR3 = vr_move1.VR3;
		VR4 = vr_move1.VR4;

		tapsearch tap_search1 = tap_searching_first(VR1, VR2, VR3, VR4, Dl_new, v0, Nnum, Z, Node_f, Lnum_a, Lnum_b, Lnum_c);
		
		mat results1;
		aVR1 = tap_search1.aVR1;
		aVR2 = tap_search1.aVR2;
		aVR3 = tap_search1.aVR3;
		aVR4 = tap_search1.aVR4;
		results1 = tap_search1.results;
		temp = results1(0, 0);
		cout << "results: " << temp << endl;
		Dl_vr(0, span(13, 15)) = aVR1;
		Dl_vr(12, span(13, 15)) = aVR2;
		Dl_vr(34, span(13, 15)) = aVR3;
		Dl_vr(153, span(13, 15)) = aVR4;
	}

	mat re1 = results(0, span::all);
	cout << "Got Results!!!!!!!!" << endl;
	//system("pause");
	vr_reg.Dl_vr = Dl_vr;
	vr_reg.re_vr = re1;
	return vr_reg;
}
