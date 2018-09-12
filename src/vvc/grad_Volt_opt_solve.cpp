# include <armadillo>
# include <iostream>
# include <math.h>
# include "fun_return.h"
# include "load_system_data.h"

using namespace std;
using namespace arma;

gradV grad_Volt_opt_solve(mat Dl_vr, double v0, cx_mat Z, double bkva, mat Dl_pv, mat PV_node)
{
	gradV grad;
	//double beta0 = 0.5;
	double beta0 = 0.5;
	double lb_v = 0.95 - 0.00001;
	double ub_v = 1.05 + 0.00001;
	double m_mx = 300;
	int trigger;
	double Ploss_afterDeVVC, Ploss_red, Qinj_endnode;
	mat Qset_a, Qset_b, Qset_c;
	mat Dl_dc;
	double Vmin_dc, Vmax_dc;
	int number1;

	sysdata sysinfo = load_system_data();
	//cx_mat Z = sysinfo.Z;

	int Ldl = Dl_vr.n_rows;
	int Wdl = Dl_vr.n_cols;

	mat Dl = Dl_vr;

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
		if ((int)Dl_pv(i, 6) != 0)
			Lla = Lla + 1;
		if ((int)Dl_pv(i, 8) != 0)
			Llb = Llb + 1;
		if ((int)Dl_pv(i, 10) != 0)
			Llc = Llc + 1;
	}
	cnt_nodes = cnt_nodes + 1;//No of nodes= No of branches +1
	mat Node_f = zeros(1, cnt_nodes);
	mat Load_a = zeros(1, Lla);
	mat Load_b = zeros(1, Llb);
	mat Load_c = zeros(1, Llc);

	cout << Lla << " " << Llb << " " << Llc << endl;
	
	for (int i = 0; i < Ldl && j<cnt_nodes && ja<Lla && jb<Llb &&jc<Llc; ++i)
	{
		if ((int)Dl(i, 0) != 0)
		{
			Node_f(0, j) = Dl(i, 2);
			++j;
		}
		/*if ((int)Dl_pv(i, 6) != 0)
		{
			Load_a(ja) = Dl_pv(i, 2);
			++ja;
		}
		if ((int)Dl_pv(i, 8) != 0)
		{
			Load_b(jb) = Dl_pv(i, 2);
			++jb;
		}
		if ((int)Dl_pv(i, 10) != 0)
		{
			Load_c(jc) = Dl_pv(i, 2);
			++jc;
		}*/
	}
	Node_f = Node_f.st();

	for (int i = 0; i < Ldl; ++i)
	{
		if ((int)Dl_pv(i, 6) != 0)
		{
			Load_a(ja) = Dl_pv(i, 2);
			++ja;
		}
		if ((int)Dl_pv(i, 8) != 0)
		{
			Load_b(jb) = Dl_pv(i, 2);
			++jb;
		}
		if ((int)Dl_pv(i, 10) != 0)
		{
			Load_c(jc) = Dl_pv(i, 2);
			++jc;
		}
	}

	//cout << Load_a.st() << endl;
	/*uvec xa, xb, xc, y;
	xa = find(Dl_pv.col(6) == PV_node);
	xb = find(Dl_pv.col(8) == PV_node);
	xc = find(Dl_pv.col(10) == PV_node);
	y << 2;
	Load_a = Dl(xa, y);
	Load_b = Dl(xb, y);
	Load_c = Dl(xc, y);
	//cout << Node_f << endl;
	cout << Load_a.st() << endl;
	cout << Load_b.st() << endl;
	cout << Load_c.st() << endl;*/
	y_re Y_return = form_Y_abc(Dl, sysinfo.Z, sysinfo.bkva, sysinfo.bkv);
	int Lbr = Y_return.brnches.n_rows;
	int Wbr = Y_return.brnches.n_cols;

	//save three phase branch data separately
	cx_mat brn_a = cx_mat(zeros(Y_return.Lnum_a, Y_return.brnches.n_cols), zeros(Y_return.Lnum_a, Y_return.brnches.n_cols));
	cx_mat brn_b = cx_mat(zeros(Y_return.Lnum_b, Y_return.brnches.n_cols), zeros(Y_return.Lnum_b, Y_return.brnches.n_cols));
	cx_mat brn_c = cx_mat(zeros(Y_return.Lnum_c, Y_return.brnches.n_cols), zeros(Y_return.Lnum_c, Y_return.brnches.n_cols));

	ja = 0;
	jb = 0;
	jc = 0;

	for (int i = 0; i < Lbr && ja<Y_return.Lnum_a && jb<Y_return.Lnum_b && jc<Y_return.Lnum_c; ++i)
	{
		if (abs(Y_return.brnches(i, 2)) != 0)
		{
			brn_a.row(ja) = Y_return.brnches.row(i);
			++ja;
		}
		if (abs(Y_return.brnches(i, 3)) != 0)
		{
			brn_b.row(jb) = Y_return.brnches.row(i);
			++jb;
		}
		if (abs(Y_return.brnches(i, 4)) != 0)
		{
			brn_c.row(jc) = Y_return.brnches.row(i);
			++jc;
		}
	}

	mat Dl_bf_vvc = Dl;

	VPQ dpf_re = DPF_return123(Dl_bf_vvc, Z);
	//cout << "Vpolar = \n" << dpf_re.Vpolar << endl;
	//cout << "PQb = \n" << dpf_re.Vpolar << endl;
	mat Vpolar = dpf_re.Vpolar;
	mat PQb = dpf_re.PQb;
	mat PQL = dpf_re.PQL;

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
	//cout << "total load (kW) per phase:" << Pload_total << endl;
	//cout << "total loss (kW):" << Ploss_orig << endl;



	//document the valid Vs and corrsponding node number for each phase
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

	//Qset before update
	Qset_a = dpf_re.Qset_a;
	Qset_b = dpf_re.Qset_b;
	Qset_c = dpf_re.Qset_c;
	mat Dl_new = Dl;
	Dl_new.col(7) = Qset_a;
	Dl_new.col(9) = Qset_b;
	Dl_new.col(11) = Qset_c;

	//newbrn Newbrn_return = rename_brn(Node_a, Node_b, Node_c, brn_a, brn_b, brn_c, Y_return.Lnum_a, Y_return.Lnum_b, Y_return.Lnum_c, Lna, Lnb, Lnc);

	//cout << "Phase A Vmin:" << min(V_a) << endl;
	//cout << "Phase B Vmin:" << min(V_b) << endl;
	//cout << "Phase C Vmin:" << min(V_c) << endl;
	mat Vmin_abc, Vmax_abc;
	Vmin_abc << min(min(V_a)) << min(min(V_b)) << min(min(V_c)) << endr;
	Vmax_abc << max(max(V_a)) << max(max(V_b)) << max(max(V_c)) << endr;
	double Vmin_orig = min(min(Vmin_abc));
	double Vmax_orig = max(max(Vmax_abc));
	//cout << "Vmax (p.u.) = " << Vmax_orig << endl;
	//cout << "Vmin (p.u.) = " << Vmin_orig << endl;

	if (Vmin_orig < lb_v || Vmax_orig > ub_v)
	{
		trigger = 1;
		
		number1 = 0;
		double Vmin, Vmax;
		mat Fv_a, Fv_b, Fv_c;
		Fv_a = zeros(Lnum_a, 1);
		Fv_b = zeros(Lnum_b, 1);
		Fv_c = zeros(Lnum_c, 1);
		Vmin = Vmin_orig;
		Vmax = Vmax_orig;

		mat V_all, histfeas;
		histfeas = zeros(m_mx, 12);
		field<mat> V_his(m_mx,1);

		while (number1 < m_mx)
		{
			if (Vmin >= lb_v && Vmax <= ub_v)
			{
				break;
			}
			else
			{
				for (int i = 1; i <= Lnum_a; ++i)
				{
					if (V_a(i) < lb_v || V_a(i) > ub_v)
					{
						Fv_a(i - 1, 0) = 2 * (V_a(i) - 1);
					}
					else
					{
						Fv_a(i - 1, 0) = 0;
					}
				}
				
				for (int i = 1; i <= Lnum_b; ++i)
				{
					if (V_b(i) < lb_v || V_b(i) > ub_v)
					{
						Fv_b(i - 1, 0) = 2 * (V_b(i) - 1);
					}
					else
					{
						Fv_b(i - 1, 0) = 0;
					}
				}

				for (int i = 1; i <= Lnum_c; ++i)
				{
					if (V_c(i) < lb_v || V_c(i) > ub_v)
					{
						Fv_c(i - 1, 0) = 2 * (V_c(i) - 1);
					}
					else
					{
						Fv_c(i - 1, 0) = 0;
					}
				}
				
				//form Jacobian Matrices for each phase
				mat L_a, L_b, L_c;
				L_a = form_L(Y_return.Y_a, V_a, theta_a, Lnum_a+1);
				L_b = form_L(Y_return.Y_b, V_b, theta_b, Lnum_b+1);
				L_c = form_L(Y_return.Y_c, V_c, theta_c, Lnum_c+1);

				mat lambda_a = -inv(L_a.st())*Fv_a;
				mat lambda_b = -inv(L_b.st())*Fv_b;
				mat lambda_c = -inv(L_c.st())*Fv_c;

				//form deltaG/deltau
				mat gu_a, gu_b, gu_c;
				gu_a = gu_a.eye(Lnum_a, Lnum_a);
				gu_b = gu_b.eye(Lnum_b, Lnum_b);
				gu_c = gu_c.eye(Lnum_c, Lnum_c);

				arma::mat g_vq_a = -gu_a.st()*lambda_a;//Gradient in p.u.  aka df/du , u is Qinj
				arma::mat g_vq_b = -gu_b.st()*lambda_b;
				arma::mat g_vq_c = -gu_c.st()*lambda_c;
				//cout << "g_vq_a" << g_vq_a << endl;
				arma::mat g_min;
				g_min << min(min(abs(g_vq_a))) << min(min(abs(g_vq_b))) << min(min(abs(g_vq_c))) << arma::endr;
				arma::mat g_max;
				g_max << max(max(abs(g_vq_a))) << max(max(abs(g_vq_b))) << max(max(abs(g_vq_c))) << arma::endr;
				double gmax = max(max(g_max));
				double gmin = min(min(g_min));
				//std::cout << "max gradient=" << gmax << "	" << "min grad=" << gmin << "(p.u.)" << std::endl;
				double gabs_min = min(min(abs(join_cols(g_vq_a, join_cols(g_vq_b, g_vq_c)))));
				// end of gradient calculation	

				mat grad_jar;
				grad_jar = max(max(abs(join_cols(g_vq_a, join_cols(g_vq_b, g_vq_c)))));
				double gabs_max = max(abs(grad_jar.elem(find(grad_jar != 0))));

				//cout << "gabs_max=" << gabs_max << endl;
				double cvq_a = beta0 / (sysinfo.bkva / 3) / gabs_max;
				double cvq_b = beta0 / (sysinfo.bkva / 3) / gabs_max;
				double cvq_c = beta0 / (sysinfo.bkva / 3) / gabs_max;
				//cout << Load_a << "\n\n" << Node_a << endl;
				double g;
				//for phase A
				for (int i = 0; i < Lla; ++i)
				{
					//cout << "Reached Here for i = " << i << endl;
					for (int j = 0; j < Lna; ++j)
					{
						//cout << "Reached Here for j = " << j << endl;
						if (Node_a(j) == Load_a(i))
						{
							g = g_vq_a(j - 1)*bkva;
							for (int m = 0; m < Ldl; ++m)
							{
								if (Dl(m, 2) == Node_a(j))
								{
									Dl(m, 7) = Dl(m, 7) + cvq_a*g;
									//cout << "m = " << m << endl;									
								}
							}
						}
					}
				}
				
				//for phase B
				for (int i = 0; i < Llb; ++i)
				{
					for (int j = 0; j < Lnb; ++j)
					{
						if (Node_b(j) == Load_b(i))
						{
							g = g_vq_b(j - 1)*bkva;
							for (int m = 0; m < Ldl; ++m)
							{
								if (Dl(m, 2) == Node_b(j))
								{
									Dl(m, 9) = Dl(m, 9) + cvq_b*g;
								}
							}
						}
					}
				}

				//for phase C
				for (int i = 0; i < Llc; ++i)
				{
					for (int j = 0; j < Lnc; ++j)
					{
						if (Node_c(j) == Load_c(i))
						{
							g = g_vq_c(j - 1)*bkva;
							for (int m = 0; m < Ldl; ++m)
							{
								if (Dl(m, 2) == Node_c(j))
								{
									Dl(m, 11) = Dl(m, 11) + cvq_c*g;
								}
							}
						}
					}
				}
				
				VPQ dpf_re = DPF_return123(Dl, Z);
				mat Vpolar = dpf_re.Vpolar;
				mat PQb = dpf_re.PQb;
				mat PQL = dpf_re.PQL;

				int Lvp = Vpolar.n_rows;
				int Wvp = Vpolar.n_cols;
				int Lpq = PQb.n_rows;
				int Wpq = PQb.n_cols;
				double Pla_t = sum(PQL.col(0));
				double Plb_t = sum(PQL.col(2));
				double Plc_t = sum(PQL.col(4));

				// calculate power loss
				mat Pload_total, Ploss_oph;
				double Ploss_osize;
				Pload_total << Pla_t << Plb_t << Plc_t << endr;
				Ploss_oph << (PQb(0, 0) - Pla_t) << (PQb(0, 2) - Plb_t) << (PQb(0, 4) - Plc_t) << endr;
				Ploss_osize = accu(Ploss_oph);
				//cout << "Ploss_osize: " << Ploss_osize << endl;

				Vabc V_abc = V_abc_list(Vpolar, Node_f, Lvp, Lnum_a, Lnum_b, Lnum_c);

				V_a = V_abc.V_a;
				V_b = V_abc.V_b;
				V_c = V_abc.V_c;

				mat Vmin_abc, Vmax_abc;
				Vmin_abc << min(min(V_a)) << min(min(V_b)) << min(min(V_c)) << endr;
				Vmax_abc << max(max(V_a)) << max(max(V_b)) << max(max(V_c)) << endr;
				Vmin = min(min(Vmin_abc));
				Vmax = max(max(Vmax_abc));

				V_all = join_cols(V_a, join_cols(V_b, V_c));
				
				V_his(number1,0) = V_all;
				mat temp;
				temp << number1 << Vmax << Vmin << Dl(31, 7) << Dl(31, 9) << Dl(31, 11) << Dl(40, 7) << Dl(40, 9) << Dl(40, 11) << gmax << gmin << cvq_a << endr;
				//histfeas(number1, span::all) = temp;
				histfeas.row(number1) = temp;
				if (Vmin > lb_v && Vmax < ub_v)
				{
					break;
				}
				number1 = number1 + 1;
			}
		}
		
		mat du = Dl(span::all, span(6, 11));
		mat du_de = du;

		//mat Dl_size = Dl_bf_vvc;
		//mat ta, tb, tc, s_sq;
		//ta = Dl_size.col(6);
		//ta.transform([](double val) { return (pow(val,2)); });
		//tb = Dl_size.col(8);
		//tb.transform([](double val) { return (pow(val, 2)); });
		//tc = Dl_size.col(10);
		//tc.transform([](double val) { return (pow(val, 2)); });
		//s_sq = join_rows(ta, join_rows(tb, tc));

		//s_sq = pow(1.25, 2)*s_sq;

		Dl_dc = Dl_bf_vvc;
		Dl_dc(span::all, span(6, 11)) = du_de;

		VPQ dpf_re = DPF_return123(Dl, Z);
		mat Vpolar = dpf_re.Vpolar;
		mat PQb = dpf_re.PQb;
		mat PQL = dpf_re.PQL;

		int Lvp = Vpolar.n_rows;
		int Wvp = Vpolar.n_cols;
		int Lpq = PQb.n_rows;
		int Wpq = PQb.n_cols;
		double Pla_t = sum(PQL.col(0));
		double Plb_t = sum(PQL.col(2));
		double Plc_t = sum(PQL.col(4));

		// calculate power loss
		mat Pload_total, Ploss_oph;
		double Ploss_osize;
		Pload_total << Pla_t << Plb_t << Plc_t << endr;
		Ploss_oph << (PQb(0, 0) - Pla_t) << (PQb(0, 2) - Plb_t) << (PQb(0, 4) - Plc_t) << endr;
		Ploss_osize = accu(Ploss_oph);

		Vabc V_abc = V_abc_list(Vpolar, Node_f, Lvp, Lnum_a, Lnum_b, Lnum_c);

		V_a = V_abc.V_a;
		V_b = V_abc.V_b;
		V_c = V_abc.V_c;

		mat Vmin_abc, Vmax_abc;
		Vmin_abc << min(min(V_a)) << min(min(V_b)) << min(min(V_c)) << endr;
		Vmax_abc << max(max(V_a)) << max(max(V_b)) << max(max(V_c)) << endr;
		Vmin = min(min(Vmin_abc));
		Vmax = max(max(Vmax_abc));

		Ploss_afterDeVVC = Ploss_osize;
		Ploss_red = Ploss_orig - Ploss_afterDeVVC;
		int l = Dl_dc.n_rows;
		Qinj_endnode = Dl_dc(l - 1, 11);

		Qset_a = Dl_dc(span::all, 7);
		Qset_b = Dl_dc(span::all, 9);
		Qset_c = Dl_dc(span::all, 11);

		//cout << "Qseta : \n" << Qset_a << endl << "Qsetb : \n" << Qset_b << endl << "Qsetc : \n" << Qset_c << endl;
		//cout << "hist_feas: " << histfeas << endl;

		Vmin_dc = Vmin;
		Vmax_dc = Vmax;
	}
	else
	{
		trigger = 0;
		number1 = 0;
		cout << "No voltage violation, do not need SSTs to regulate voltage" << endl;
		Dl_dc = Dl_bf_vvc;
		Vmin_dc = Vmin_orig;
		Vmax_dc = Vmax_orig;
		Ploss_afterDeVVC = Ploss_orig;
	}

	grad.Dl_volt = Dl_dc;
	grad.Ploss_afterVR = Ploss_afterDeVVC;
	grad.Vmax = Vmax_dc;
	grad.Vmin = Vmin_dc;
	grad.trigger1 = trigger;
	grad.no_itea = number1;
	return grad;
}
