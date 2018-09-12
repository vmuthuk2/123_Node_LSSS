# include <armadillo>
# include <iostream>
# include <math.h>
# include "fun_return.h"
# include "load_system_data.h"

using namespace std;
using namespace arma;

gradVAR var_opt_loss_solver(arma::mat Dl_vr, double v0, arma::cx_mat Z, double bkva, arma::mat Dl_pv, arma::mat PV_node)
{
	gradVAR grad;
	int itea_stepsch = 0;

	double beta0 = 0.02;
	double alpha = 1.1;
	double lb_v = 0.95;
	double ub_v = 1.05;
	double m_max = 300;
	double epsi = 0.0001;
	
	mat Dl = Dl_vr;
	int Ldl = Dl.n_rows;
	int Wdl = Dl.n_cols;

	//int Ldl = Dl_vr.n_rows;
	//int Wdl = Dl_vr.n_cols;

	//mat Dl = Dl_vr;
	sysdata sysinfo = load_system_data();
	
	mat du, step_size; // delta control
	double Ploss_aftter_ctrl, Ploss_red;
	bool flag = true;
	double Vmax, Vmin;

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

	for (int i = 0; i < Ldl && j<cnt_nodes && ja<Lla && jb<Llb && jc<Llc; ++i)
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
	/*uvec xa, xb, xc, y;
	xa = find(Dl.col(6) != 0);
	xb = find(Dl.col(8) != 0);
	xc = find(Dl.col(10) != 0);
	y << 2;
	Load_a = Dl(xa, y);
	Load_b = Dl(xb, y);
	Load_c = Dl(xc, y);
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
	mat Qset_a = dpf_re.Qset_a;
	mat Qset_b = dpf_re.Qset_b;
	mat Qset_c = dpf_re.Qset_c;
	mat Dl_new = Dl;

	Dl_new.col(7) = Qset_a;
	Dl_new.col(9) = Qset_b;
	Dl_new.col(11) = Qset_c;

	newbrn Newbrn_return = rename_brn(Node_a, Node_b, Node_c, brn_a, brn_b, brn_c, Y_return.Lnum_a, Y_return.Lnum_b, Y_return.Lnum_c, Lna, Lnb, Lnc);

	//cout << "Phase A Vmin:" << min(V_a) << endl;
	//cout << "Phase B Vmin:" << min(V_b) << endl;
	//cout << "Phase C Vmin:" << min(V_c) << endl;
	mat Vmin_abc, Vmax_abc;
	Vmin_abc << min(min(V_a)) << min(min(V_b)) << min(min(V_c)) << endr;
	Vmax_abc << max(max(V_a)) << max(max(V_b)) << max(max(V_c)) << endr;
	double Vmin_orig = min(min(Vmin_abc));
	double Vmax_orig = max(max(Vmax_abc));
	
	//deltaF/deltaTheta
	arma::mat Ftheta_a = form_Ftheta(Y_return.Y_a, V_a, theta_a, Newbrn_return.newbrn_a, Lna, Lnum_a);
	arma::mat Ftheta_b = form_Ftheta(Y_return.Y_b, V_b, theta_b, Newbrn_return.newbrn_b, Lnb, Lnum_b);
	arma::mat Ftheta_c = form_Ftheta(Y_return.Y_c, V_c, theta_c, Newbrn_return.newbrn_c, Lnc, Lnum_c);

	//detltaF/deltaV
	arma::mat Fv_a = form_Fv(Y_return.Y_a, V_a, theta_a, Newbrn_return.newbrn_a, Lna, Lnum_a);
	arma::mat Fv_b = form_Fv(Y_return.Y_b, V_b, theta_b, Newbrn_return.newbrn_b, Lnb, Lnum_b);
	arma::mat Fv_c = form_Fv(Y_return.Y_c, V_c, theta_c, Newbrn_return.newbrn_c, Lnc, Lnum_c);
	//detltaF/deltaX
	arma::mat Fx_a = join_cols(Ftheta_a, Fv_a);
	arma::mat Fx_b = join_cols(Ftheta_b, Fv_b);
	arma::mat Fx_c = join_cols(Ftheta_c, Fv_c);


	//form Jacobian Matrices for each phase
	arma::mat J_a, J_b, J_c;
	J_a = form_J(Y_return.Y_a, V_a, theta_a, Lna);
	J_b = form_J(Y_return.Y_b, V_b, theta_b, Lnb);
	J_c = form_J(Y_return.Y_c, V_c, theta_c, Lnc);

	// get lambda for each phase
	arma::mat lambda_a = -inv(J_a.st())*Fx_a;//need LAPAC etc
	arma::mat lambda_b = -inv(J_b.st())*Fx_b;
	arma::mat lambda_c = -inv(J_c.st())*Fx_c;

	//std::cout << lambda_a << std::endl;

	//deltaG/delta_Qnj from sst
	//deltaP/delta_Qinj
	arma::mat Gpq_a = arma::zeros(Lnum_a, Lla);
	arma::mat Gpq_b = arma::zeros(Lnum_b, Llb);
	arma::mat Gpq_c = arma::zeros(Lnum_c, Llc);
	//deltaQ/delta_Qinj
	arma::mat Gqq_a = arma::zeros(Lnum_a, Lla);
	arma::mat Gqq_b = arma::zeros(Lnum_b, Llb);
	arma::mat Gqq_c = arma::zeros(Lnum_c, Llc);
	
	//Load_a = Load_a.st();
	//Load_b = Load_b.st();
	//Load_c = Load_c.st();
	//Phase A
	for (ia = 0; ia < Lnum_a; ++ia)
	{
		for (ja = 0; ja < Lla; ++ja)
		{
			if (Node_a(0, ia + 1) == Load_a(0, ja))
			{
				Gqq_a(ia, ja) = -1;
			}

		}
	}

	//Phase B
	for (ib = 0; ib < Lnum_b; ++ib)
	{
		for (jb = 0; jb < Llb; ++jb)
		{
			if (Node_b(0, ib + 1) == Load_b(0, jb))
			{
				Gqq_b(ib, jb) = -1;
			}

		}
	}

	//Phase C
	for (ic = 0; ic < Lnum_c; ++ic)
	{
		for (jc = 0; jc < Llc; ++jc)
		{
			if (Node_c(0, ic + 1) == Load_c(0, jc))
			{
				Gqq_c(ic, jc) = -1;
			}

		}
	}

	//form deltaG/deltaQinj
	arma::mat gu_a, gu_b, gu_c;
	gu_a = join_cols(Gpq_a, Gqq_a);
	gu_b = join_cols(Gpq_b, Gqq_b);
	gu_c = join_cols(Gpq_c, Gqq_c);



	arma::mat g_vq_a = -gu_a.st()*lambda_a;//Gradient in p.u.  aka df/du , u is Qinj
	arma::mat g_vq_b = -gu_b.st()*lambda_b;
	arma::mat g_vq_c = -gu_c.st()*lambda_c;
	//std::cout<< g_vq_a << std::endl;

	arma::mat g_min;
	g_min << min(min(abs(g_vq_a))) << min(min(abs(g_vq_b))) << min(min(abs(g_vq_c))) << arma::endr;
	arma::mat g_max;
	g_max << max(max(abs(g_vq_a))) << max(max(abs(g_vq_b))) << max(max(abs(g_vq_c))) << arma::endr;
	double gmax = max(max(g_max));
	double gmin = min(min(g_min));
	//std::cout << "max gradient=" << gmax << "	" << "min grad=" << gmin << "(p.u.)" << std::endl;
	double gabs_min = min(min(abs(join_cols(g_vq_a, join_cols(g_vq_b, g_vq_c)))));
	// end of gradient calculation	

	//cout << "gabs_min=" << gabs_min << endl;
	double cvq_a = beta0 / (sysinfo.bkva / 3) / gabs_min;
	double cvq_b = beta0 / (sysinfo.bkva / 3) / gabs_min;
	double cvq_c = beta0 / (sysinfo.bkva / 3) / gabs_min;

	mat ctrl_o = Dl;// save the previous control


	for (int m = 0; m < m_max; m++)
	{
		// update Qinj for three phase separately
		double Gupdate;
		//Phase A
		for (ia = 0; ia < Lla; ++ia)
		{
			for (ja = 0; ja < Ldl; ++ja)
			{
				if (Dl(ja, 2) == Load_a(0, ia))
				{
					Gupdate = g_vq_a(ia, 0)*(sysinfo.bkva / 3)*cvq_a;
					//Dl_new(ja, 7) = ctrl_o(ja, 7) - Gupdate;
					Dl(ja, 7) = ctrl_o(ja, 7) - Gupdate;
				}
			}
		}//end of Phase A
		//cout << "Gupdate: " << Dl_new.col(7) << endl;

		 //Phase B
		for (ib = 0; ib < Llb; ++ib)
		{
			for (jb = 0; jb < Ldl; ++jb)
			{
				if (Dl(jb, 2) == Load_b(0, ib))
				{
					Gupdate = g_vq_b(ib, 0)*(sysinfo.bkva / 3)*cvq_b;
					//Dl_new(jb, 9) = ctrl_o(jb, 9) - Gupdate;
					Dl(jb, 9) = ctrl_o(jb, 9) - Gupdate;
				}
			}
		}//end of Phase B

		 //Phase C
		for (ic = 0; ic < Llc; ++ic)
		{
			for (jc = 0; jc < Ldl; ++jc)
			{
				if (Dl(jc, 2) == Load_c(0, ic))
				{
					Gupdate = g_vq_c(ic, 0)*(sysinfo.bkva / 3)*cvq_c;
					//Dl_new(jc, 11) = ctrl_o(jc, 11) - Gupdate;
					Dl(jc, 11) = ctrl_o(jc, 11) - Gupdate;
				}
			}
		}// end of Phase C
		//mat Dl_osize = Dl_new;
		mat Dl_osize = Dl;
		mat du_temp = Dl_osize - ctrl_o;
		du = du_temp.cols(6, 11);

		// cout << "Dl_new = \n" << Dl_new << endl;
		// DPF based on Dl_new

		VPQ dpf_re = DPF_return123(Dl_osize, Z);
		//cout << "Vpolar = \n" << dpf_re.Vpolar << endl;
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
		//cout << "total load (kW) per phase:" << Pload_total << endl;
		//cout << "total loss (kW):" << Ploss_osize << endl;

		Vabc V_abc = V_abc_list(Vpolar, Node_f, Lvp, Lnum_a, Lnum_b, Lnum_c);


		V_a = V_abc.V_a;
		V_b = V_abc.V_b;
		V_c = V_abc.V_c;

		mat Vmin_abc, Vmax_abc;
		Vmin_abc << min(min(V_a)) << min(min(V_b)) << min(min(V_c)) << endr;
		Vmax_abc << max(max(V_a)) << max(max(V_b)) << max(max(V_c)) << endr;
		Vmin = min(min(Vmin_abc));
		Vmax = max(max(Vmax_abc));


		step_size << cvq_a << cvq_b << cvq_c << endr;
		//cout << "\n \n step-size at the " << m + 1 << "th iteration = " << step_size(0, 0) << endl;

		// new step-size
		cvq_a = alpha*cvq_a;
		cvq_b = alpha*cvq_b;
		cvq_c = alpha*cvq_c;

		//update Dl_new based on the new step-size
		//Phase A
		for (ia = 0; ia < Lla; ++ia)
		{
			for (ja = 0; ja < Ldl; ++ja)
			{
				if (Dl(ja, 2) == Load_a(0, ia))
				{
					Gupdate = g_vq_a(ia, 0)*(sysinfo.bkva / 3)*cvq_a;
					//Dl_new(ja, 7) = ctrl_o(ja, 7) - Gupdate;
					Dl(ja, 7) = ctrl_o(ja, 7) - Gupdate;
				}
			}
		}//end of Phase A

		 //Phase B
		for (ib = 0; ib < Llb; ++ib)
		{
			for (jb = 0; jb < Ldl; ++jb)
			{
				if (Dl(jb, 2) == Load_b(0, ib))
				{
					Gupdate = g_vq_b(ib, 0)*(sysinfo.bkva / 3)*cvq_b;
					//Dl_new(jb, 9) = ctrl_o(jb, 9) - Gupdate;
					Dl(jb, 9) = ctrl_o(jb, 9) - Gupdate;
				}
			}
		}//end of Phase B

		 //Phase C
		for (ic = 0; ic < Llc; ++ic)
		{
			for (jc = 0; jc < Ldl; ++jc)
			{
				if (Dl(jc, 2) == Load_c(0, ic))
				{
					Gupdate = g_vq_c(ic, 0)*(sysinfo.bkva / 3)*cvq_c;
					//Dl_new(jc, 11) = ctrl_o(jc, 11) - Gupdate;
					Dl(jc, 11) = ctrl_o(jc, 11) - Gupdate;
				}
			}
		}// end of Phase C

		//mat Dl_nsize = Dl_new;
		mat Dl_nsize = Dl;

		dpf_re = DPF_return123(Dl_nsize, Z);// run DPF based on new step-size
		Vpolar = dpf_re.Vpolar;
		PQb = dpf_re.PQb;
		PQL = dpf_re.PQL;
		//cout << "Vpolar: \n" << Vpolar << endl;
		//system("pause");

		Pla_t = sum(PQL.col(0));
		Plb_t = sum(PQL.col(2));
		Plc_t = sum(PQL.col(4));

		// calculate power loss
		mat Ploss_nph;
		double Ploss_nsize;
		Pload_total << Pla_t << Plb_t << Plc_t << endr;
		Ploss_nph << (PQb(0, 0) - Pla_t) << (PQb(0, 2) - Plb_t) << (PQb(0, 4) - Plc_t) << endr;
		Ploss_nsize = accu(Ploss_nph);
		//cout << "total new loss (kW):" << Ploss_nsize << endl;
		//cout << "total old loss (kW):" << Ploss_osize << endl;

		if (Ploss_nsize>Ploss_osize)
		{
			Dl = Dl_osize;//save the Dl based on old size
			Ploss_aftter_ctrl = Ploss_osize;
			itea_stepsch = m+1;
			cout << " Best step-size obtained! Searching ends at the " << m + 1 << "th iteration \n";
			cout << " Expected power loss (kW) = " << Ploss_aftter_ctrl << endl;
			cout << " Expected power loss reduction (kW) = " << Ploss_orig - Ploss_aftter_ctrl << endl;
			cout << " Actual Step_Size: " << step_size << endl;
			//step_size;
			break;
		}
		else
		{
			Ploss_aftter_ctrl = Ploss_osize;
			itea_stepsch = m+1;
		}

		if (m == m_max)
		{
			//cout << " Unable to obtain the best step-size withing " << m + 1 << "th iterations \n";
			itea_stepsch = m;
			break;
		}

		if ((Vmax > 1.05 + epsi) || (Vmin < 0.95 - epsi))
		{
			cout << "Exit Gradient based Ploss Minimization due to voltage violation!" << endl;
			Ploss_aftter_ctrl = Ploss_orig;
			//Ploss_red = 0;
			Vmin = Vmin_orig;
			Vmax = Vmax_orig;
			Qset_a = Dl_bf_vvc.col(7);
			Qset_b = Dl_bf_vvc.col(9);
			Qset_c = Dl_bf_vvc.col(11);
			cout << " Best step-size obtained! Searching ends at the " << m + 1 << "th iteration \n";
			cout << " Expected power loss (kW) = " << Ploss_aftter_ctrl << endl;
			cout << " Expected power loss reduction (kW) = " << Ploss_orig - Ploss_aftter_ctrl << endl;
			cout << " Actual Step_Size: " << step_size << endl;
			break;
		}
	}// end of step-size search

	if (Ploss_aftter_ctrl > Ploss_orig)
	{
		cout << " Direction of gradients need to be reversed! \n";
		flag = false;
		Vmin = Vmin_orig;
		Vmax = Vmax_orig;
		Ploss_aftter_ctrl = Ploss_orig;
		Ploss_red = 0;
		Qset_a = Dl_bf_vvc.col(7);
		Qset_b = Dl_bf_vvc.col(9);
		Qset_c = Dl_bf_vvc.col(11);
	}
	else
	{
		mat du_de = du;
		mat Dl_dc = Dl_bf_vvc;
		Dl_dc(span::all, span(6, 11)) = Dl_dc(span::all, span(6, 11)) + du_de;

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
		//double Ploss_orig;
		Pload_total << Pla_t << Plb_t << Plc_t << endr;
		Ploss_orig_ph << (PQb(0, 0) - Pla_t) << (PQb(0, 2) - Plb_t) << (PQb(0, 4) - Plc_t) << endr;
		Ploss_aftter_ctrl = accu(Ploss_orig_ph);

		Ploss_red = Ploss_aftter_ctrl - Ploss_orig;
		//Qset_a = dpf_re.Qset_a;
		//Qset_b = dpf_re.Qset_b;
		//Qset_c = dpf_re.Qset_c;
		Qset_a = Dl_dc.col(7);
		Qset_b = Dl_dc.col(9);
		Qset_c = Dl_dc.col(11);
	}
	//system("pause");
	grad.V_mn = Vmin;
	grad.V_mx = Vmax;
	//grad.Qcmd = join_cols(Qset_a, join_cols(Qset_b, Qset_c));
	grad.Qcmd = join_rows(Qset_a, join_rows(Qset_b, Qset_c));
	grad.totalloss = Ploss_aftter_ctrl;
	grad.itea_stepsch = itea_stepsch;
	return grad;
}
