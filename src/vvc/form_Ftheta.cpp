//# define ARMA_DONT_USE_WRAPPER

#include "armadillo"
#include "fun_return.h"
#include "cmath"
#include <iostream>

using namespace arma;

arma::mat form_Ftheta(arma::cx_mat Y, arma::mat V, arma::mat theta, arma::cx_mat brn, int Ln, int Lnm)
{
	//cout << brn.col(0) << " " << brn.col(1) << endl;
	//cout << "Y: " << Y.n_rows << "*" << Y.n_cols << endl;
	//cout << "V: " << V.n_rows << "*" << V.n_cols << endl;
	//cout << "theta: " << theta.n_rows << "*" << theta.n_cols << endl;
	//cout << "brn: " << brn.n_rows << "*" << endl;
	//cout << "Ln: " << Ln << endl;
	//cout << "Lnm: " << Lnm << endl;
	//cout << "Y matrix: " << Y << endl;
	arma::mat Ftheta = arma::zeros(Ln - 1,1);
	const double pi = 3.14159265358979323846;
	int i = 0;
	int j = 0;
	long double R;
	for (i = 0; i < (Ln - 1); ++i)
	{
		R = 0;
		for (j = 0; j < Lnm; ++j)
		{//cout << "j=:"<< j << endl;
			int s = (int)real(brn(j, 0));
			int r = (int)real(brn(j, 1));
			//cout << "s=:" << s << "	" << "r=" << r << "i+1" << i+1 << endl;
			if (s == i+1)
			{
				//cout << "Y(s,r)=" << Y(s, r) << endl;
				//cout << "V(s)=" << V(s, 0) << "	" << "V(r)=" << V(r, 0) << endl;
				//cout << "theta(s)=" << theta(0, s) << "	" << "theta(r)=" << theta(0, r) << endl;
				R = R - 2 * (-real(Y(s, r)))*V(s,0)*V(r,0)*(-sin((theta(0,s) - theta(0,r))*pi/ 180));
			}
			if (r == i+1)
			{
				//cout << "Y(s,r)=" << Y(s, r) << endl;
				//cout << "V(s)=" << V(s, 0) << "	" << "V(r)=" << V(r, 0) << endl;
				//cout << "theta(s)=" << theta(0, s) << "	" << "theta(r)=" << theta(0, r) << endl;
				R = R - 2 * (-real(Y(s, r)))*V(s,0)*V(r,0)*sin((theta(0,s) - theta(0,r))*pi / 180);
			}
			//cout << "R=" << R << endl;
		}
		Ftheta(i,0) = R;
		//cout << "Ftheta(" << i << ")=" << Ftheta(i, 0) << endl;
	}
	return Ftheta;
}