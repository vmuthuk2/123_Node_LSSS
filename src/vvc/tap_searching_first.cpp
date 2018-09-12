# include <armadillo>
# include <iostream>
# include <math.h>
# include "fun_return.h"
#include "load_system_data.h"

using namespace std;
using namespace arma;

tapsearch tap_searching_first(mat VR1, mat VR2, mat VR3, mat VR4, mat Dl_bf_voltc, double v0, int Nnum, cx_mat Z, mat Node_f, int Lnum_a, int Lnum_b, int Lnum_c)
{
	cout << "Entering tap_searching_first.cpp" << endl;
	tapsearch tap_search;
	double epsi = 0.0001;
	field<mat> VRlist_(4,1);
	mat ls_len(4,1);

	mat lcd_vr1, lcd_vr2, lcd_vr3, lcd_vr4;
	lcd_vr1 << 1 << 1 << 1 << endr;
	lcd_vr2 << 1 << 1 << 1 << endr;
	lcd_vr3 << 1 << 0 << 0 << endr;
	lcd_vr4 << 1 << 0 << 1 << endr;

	field<mat> vr_lcd;
	vr_lcd << lcd_vr1 << lcd_vr2 << lcd_vr3 << lcd_vr4;

	mat a1, a2, a3, a4;
	field<mat> alpha;
	a1 = Dl_bf_voltc(0, span(13, 15));
	a2 = Dl_bf_voltc(12, span(13, 15));
	a3 = Dl_bf_voltc(34, span(13, 15));
	a4 = Dl_bf_voltc(153, span(13, 15));

	alpha << a1 << a2 << a3 << a4 << endr;
	//cout << "Alpha: " << alpha << endl;
	mat results, re1, re2;
	mat aVR1, aVR2, aVR3, aVR4;
	mat feasible_taps;

	if (accu(VR1) == 0 && accu(VR2) == 0 && accu(VR3) && accu(VR4))
	{
		results = zeros<mat>(1,16);
		results(0, span(4, 6)) = Dl_bf_voltc(0, span(13, 15));
		results(0, span(7, 9)) = Dl_bf_voltc(12, span(13, 15));
		results(0, span(10, 12)) = Dl_bf_voltc(34, span(13, 15));
		results(0, span(13, 15)) = Dl_bf_voltc(153, span(13, 15));
		cout << "Voltage Acceptable. No Search Needed" << endl;
		VPQ dpf_re = DPF_return123(Dl_bf_voltc, Z);
		mat Vpolar = dpf_re.Vpolar;
		int Lvp = Vpolar.n_rows;

		Vabc V_abc = V_abc_list(Vpolar, Node_f, Lvp, Lnum_a, Lnum_b, Lnum_c);


		mat V_a = V_abc.V_a;
		mat V_b = V_abc.V_b;
		mat V_c = V_abc.V_c;

		mat Vmin_abc, Vmax_abc;
		Vmin_abc << min(min(V_a)) << min(min(V_b)) << min(min(V_c)) << endr;
		Vmax_abc << max(max(V_a)) << max(max(V_b)) << max(max(V_c)) << endr;
		double Vmin = min(min(Vmin_abc));
		double Vmax = max(max(Vmax_abc));
		mat Vresult;
		Vresult << Vmin << Vmax << endr;
		results(0, span(2, 3)) = Vresult;
	}
	else
	{
		double alb, aub;
		alb = 0.9 - epsi;
		aub = 1.1 + epsi;

		field<mat> VR;
		VR << VR1 << VR2 << VR3 << VR4 << endr;

		mat ContainerVR, t1, t2, t3, t4, t5;

		for (int vr_idx = 0; vr_idx < VR.n_cols; ++vr_idx)
		{
			uvec temp1 = find(alpha(vr_idx) < alb, 1);
			uvec temp2 = find(alpha(vr_idx) > aub, 1);

			if (!temp1.is_empty())
			{
				cout << "Error: VR " << vr_idx+1 << " initial tap position is too low!" << endl;
				cout << alpha( vr_idx ) << endl;
				cout << "Please fix it before searching for feasible tap positions" << endl;
				break;
			}

			if (!temp2.is_empty())
			{
				cout << "Error: VR " << vr_idx+1 << " initial tap position is too high!" << endl;
				cout << alpha(vr_idx) << endl;
				cout << "Please fix it before searching for feasible tap positions" << endl;
				break;
			}

			if (temp1.is_empty() && temp2.is_empty())
			{
				cout << "VR " << vr_idx+1 << " initial tap position is feasible. Continued.... " << endl;
				t5 = alpha(vr_idx);
				cout << "size of t5: " << t5.n_rows << " * " << t5.n_cols << endl;
				//t1 = t5.transform([](double val) {return (val + (0.1 / 16)); });
				t1 = t5 + (0.1 / 16);
				t5 = alpha(vr_idx);
				//t2 = t5.transform([](double val) {return (val - (0.1 / 16)); });
				t2 = t5 - (0.1 / 16);
				t5 = alpha(vr_idx);
				//t3 = t5.transform([](double val) {return (val + (0.1 / 16)*2); });
				t3 = t5 + (0.1 / 16) * 2;
				t5 = alpha(vr_idx);
				//t4 = t5.transform([](double val) {return (val - (0.1 / 16)*2); });
				t4 = t5 - (0.1 / 16) * 2;
				
				ContainerVR = join_cols(alpha(vr_idx), t1);
				ContainerVR = join_cols(ContainerVR, t2);
				ContainerVR = join_cols(ContainerVR, t3);
				ContainerVR = join_cols(ContainerVR, t4);
			}
			
			mat VR_ratio = alpha(vr_idx);
			int rnum = ContainerVR.n_rows;
			int cnum = ContainerVR.n_cols;
			mat vr_ph = vr_lcd(vr_idx);
			mat VR_adj = VR(vr_idx);
			for (int i = 0; i < rnum; ++i)
			{
				for (int j = 0; j < cnum; ++j)
				{
					if (ContainerVR(i, j) != VR_ratio(j))
					{
						ContainerVR(i, j) = vr_ph(j)*ContainerVR(i, j);
					}

					if (ContainerVR(i, j) > aub)
					{
						ContainerVR(i, j) = 0;
					}
					else if (ContainerVR(i, j) < alb)
					{
						ContainerVR(i, j) = 0;
					}
					else
					{ }

					if (VR_adj(j) == 1 && ContainerVR(i, j) < VR_ratio(j))
					{
						ContainerVR(i, j) = 0;
					}
					if(VR_adj(j) == -1 && ContainerVR(i, j) > VR_ratio(j))
					{
						ContainerVR(i, j) = 0;
					}
					if (VR_adj(j) == 0 && ContainerVR(i, j) != VR_ratio(j))
					{
						ContainerVR(i, j) = 0;
					}
				}
			}
			//cout << "ContainerVR: \n" << ContainerVR << endl;
			//mat ContVRa = ContainerVR.col(0);
			//mat ContVRb = ContainerVR.col(1);
			//mat ContVRc = ContainerVR.col(2);

			//for (int i = 0; i < rnum; ++i)
			//{
			//	if (ContVRa(i) == 0)
			//	{
			//		ContVRa(i) = NAN;
			//	}
			//	if (ContVRb(i) == 0)
			//	{
			//		ContVRb(i) = NAN;
			//	}
			//	if (ContVRc(i) == 0)
			//	{
			//		ContVRc(i) = NAN;
			//	}
			//}

			uvec idxa, idxb, idxc;
			idxa = find(ContainerVR.col(0) != 0);
			idxb = find(ContainerVR.col(1) != 0);
			idxc = find(ContainerVR.col(2) != 0);
			mat ContVRa(idxa.n_rows,1) , ContVRb(idxb.n_rows, 1), ContVRc(idxc.n_rows, 1);

			int ia, ib, ic;
			ia = 0; ib = 0; ic = 0;
			for (int i = 0; i < rnum; ++i)
			{
				if (ContainerVR(i, 0) != 0)
				{
					ContVRa(ia, 0) = ContainerVR(i, 0);
					ia = ia + 1;
				}
				if (ContainerVR(i, 1) != 0)
				{
					ContVRb(ib, 0) = ContainerVR(i, 1);
					ib = ib + 1;
				}
				if (ContainerVR(i, 2) != 0)
				{
					ContVRc(ic, 0) = ContainerVR(i, 2);
					ic = ic + 1;
				}
			}
			//cout << "ContVRa: \n" << ContVRa << endl;
			int m = 0;
			double length_vrlist = ContVRa.n_rows * ContVRb.n_rows * ContVRc.n_rows;
			mat VR_list = zeros<mat>(length_vrlist, 3);
			for (int i = 0; i < ContVRa.n_rows; ++i)
			{
				for (int j = 0; j < ContVRb.n_rows; ++j)
				{
					for (int k = 0; k < ContVRc.n_rows; ++k)
					{
						mat t6;
						t6 << ContVRa(i) << ContVRb(j) << ContVRc(k) << endr;
						VR_list(m, span::all) = t6;
						m = m + 1;
					}
				}
			}

			VRlist_(vr_idx,0) = VR_list;
			ls_len(vr_idx) = size(VR_list, 0);
		}
		//cout << "VRlist_: \n" << VRlist_ << endl;
		mat Dl = Dl_bf_voltc;
		mat VR_ls4, VR_ls3, VR_ls2, VR_ls1;
		VR_ls4 = VRlist_(3,0);
		VR_ls3 = VRlist_(2,0);
		VR_ls2 = VRlist_(1,0);
		VR_ls1 = VRlist_(0,0);
		int cnt = 0;
		//cout << "VR_ls2: \n" << VR_ls2 << endl;
		
		for (int p = 0; p < ls_len(1); ++p)
		{
			Dl(12, span(13, 15)) = VR_ls2(p, span::all);
			for (int m = 0; m < ls_len(3); ++m)
			{
				Dl(153, span(13, 15)) = VR_ls4(m, span::all);
				for (int n = 0; n < ls_len(2); ++n)
				{
					Dl(34, span(13, 15)) = VR_ls3(n, span::all);
					for (int q = 0; q < ls_len(0); ++q)
					{
						Dl(0, span(13, 15)) = VR_ls1(q, span::all);
						VPQ dpf_temp = DPF_return123(Dl, Z);
						mat Vpolar = dpf_temp.Vpolar;
						int Lvp = Vpolar.n_rows;
						//cout << "VR1: " << VR_ls1(q, span::all) << endl;
						//cout << "VR2: " << VR_ls3(n, span::all) << endl;
						//cout << "VR3: " << VR_ls4(m, span::all) << endl;
						//cout << "VR4: " << VR_ls2(p, span::all) << endl;
						//cout << "Vpolar: \n" << Vpolar << endl;
						//system("pause");
						Vabc V_abc = V_abc_list(Vpolar, Node_f, Lvp, Lnum_a, Lnum_b, Lnum_c);

						mat V_a = V_abc.V_a;
						mat V_b = V_abc.V_b;
						mat V_c = V_abc.V_c;

						mat Vmin_abc, Vmax_abc;
						Vmin_abc << min(min(V_a)) << min(min(V_b)) << min(min(V_c)) << endr;
						Vmax_abc << max(max(V_a)) << max(max(V_b)) << max(max(V_c)) << endr;
						cout << "Phase A minimum voltage in p.u. is " << min(min(V_a)) << endl;
						cout << "Phase B minimum voltage in p.u. is " << min(min(V_b)) << endl;
						cout << "Phase C minimum voltage in p.u. is " << min(min(V_c)) << endl;
						cout << "Phase A maximum voltage in p.u. is " << max(max(V_a)) << endl;
						cout << "Phase B maximum voltage in p.u. is " << max(max(V_b)) << endl;
						cout << "Phase C maximum voltage in p.u. is " << max(max(V_c)) << endl;
						double Vmin = min(min(Vmin_abc));
						double Vmax = max(max(Vmax_abc));
						
						mat Vresult, Vpolar_a, Vpolar_b, Vpolar_c, Vpolar_temp;
						Vresult << Vmin << Vmax << endr;
						
						Vpolar_temp = Vpolar.col(0);
						Vpolar_a = Vpolar_temp.elem(find(Vpolar.col(0) != 0));
						Vpolar_temp = Vpolar.col(2);
						Vpolar_b = Vpolar_temp.elem(find(Vpolar.col(2) != 0));
						Vpolar_temp = Vpolar.col(4);
						Vpolar_c = Vpolar_temp.elem(find(Vpolar.col(4) != 0));
						
						mat VV, V_var_temp;
						VV = join_cols(Vpolar_a, Vpolar_b);
						VV = join_cols(VV, Vpolar_c);
						V_var_temp = var(VV);
						double V_var = V_var_temp(0, 0);

						//cnt = cnt + 1;
						if (Vmin >= (0.95 - epsi) && Vmax <= (1.05 + epsi))
						{
							mat temp;
							temp << V_var << cnt + 1 << Vmin << Vmax << endr;
							feasible_taps.resize(cnt + 1, 16);
							feasible_taps(cnt, span(0, 3)) = temp;
							feasible_taps(cnt, span(4, 6)) = VR_ls1.row(q);
							feasible_taps(cnt, span(7, 9)) = VR_ls2.row(p);
							feasible_taps(cnt, span(10, 12)) = VR_ls3.row(n);
							feasible_taps(cnt, span(13, 15)) = VR_ls4.row(m);
						}
						else if (Vmin < (0.95 - epsi))
						{
							mat temp;
							temp << 2*pow(10,6) << cnt + 1 << Vmin << Vmax << endr;
							feasible_taps.resize(cnt + 1, 16);
							feasible_taps(cnt, span(0, 3)) = temp;
							feasible_taps(cnt, span(4, 6)) = VR_ls1.row(q);
							feasible_taps(cnt, span(7, 9)) = VR_ls2.row(p);
							feasible_taps(cnt, span(10, 12)) = VR_ls3.row(n);
							feasible_taps(cnt, span(13, 15)) = VR_ls4.row(m);
						}
						else if (Vmax > (1.05 + epsi))
						{
							mat temp;
							temp << 1 * pow(10, 6) << cnt + 1 << Vmin << Vmax << endr;
							feasible_taps.resize(cnt + 1, 16);
							feasible_taps(cnt, span(0, 3)) = temp;
							feasible_taps(cnt, span(4, 6)) = VR_ls1.row(q);
							feasible_taps(cnt, span(7, 9)) = VR_ls2.row(p);
							feasible_taps(cnt, span(10, 12)) = VR_ls3.row(n);
							feasible_taps(cnt, span(13, 15)) = VR_ls4.row(m);
						}
						else
						{ }

						cnt = cnt + 1;

					}
				}
			}
		}

		uvec t = sort_index(feasible_taps.col(0));
		uvec cols = linspace<uvec>(0, feasible_taps.n_cols - 1, feasible_taps.n_cols);
		results = feasible_taps(t, cols);
		//cout << "Results: \n" << results << endl;
	}

	aVR1 = results(0, span(4, 6));
	aVR2 = results(0, span(7, 9));
	aVR3 = results(0, span(10, 12));
	aVR4 = results(0, span(13, 15));

	if (results(0, 0) == 1 * pow(10, 6))
	{
		uvec t = sort_index(feasible_taps.col(3));
		uvec cols = linspace<uvec>(0, feasible_taps.n_cols - 1, feasible_taps.n_cols);
		re1 = feasible_taps(t, cols);
		aVR1 = re1(0, span(4, 6));
		aVR2 = re1(0, span(7, 9));
		aVR3 = re1(0, span(10, 12));
		aVR4 = re1(0, span(13, 15));
	}
	else if (results(0, 0) == 2 * pow(10, 6))
	{
		uvec t = sort_index(feasible_taps.col(2),"descend");
		uvec cols = linspace<uvec>(0, feasible_taps.n_cols - 1, feasible_taps.n_cols);
		re2 = feasible_taps(t, cols);
		//cout << "re2: \n" << re2(span::all,span(2,6)) << endl;
		aVR1 = re2(0, span(4, 6));
		aVR2 = re2(0, span(7, 9));
		aVR3 = re2(0, span(10, 12));
		aVR4 = re2(0, span(13, 15));
	}

	cout << "Exiting tap_searching_first.cpp" << endl;
	//cout << "feasible_taps: " << feasible_taps(span::all, span(2,9)) << endl;
	//cout << "aVR1: \n" << aVR1 << endl;
	//cout << "aVR2: \n" << aVR2 << endl;
	//cout << "aVR3: \n" << aVR3 << endl;
	//cout << "aVR4: \n" << aVR4 << endl;
	//cout << "Feasible_taps: \n" << feasible_taps.col(0);
	//system("pause");
	tap_search.aVR1 = aVR1;
	tap_search.aVR2 = aVR2;
	tap_search.aVR3 = aVR3;
	tap_search.aVR4 = aVR4;
	tap_search.results = results;
	return tap_search;
}