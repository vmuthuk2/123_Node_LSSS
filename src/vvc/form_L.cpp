//# define ARMA_DONT_USE_WRAPPER

#include "armadillo"
#include "fun_return.h"
#include "cmath"
#include <iostream>

//using namespace arma;

arma::mat form_L(arma::cx_mat Y, arma::mat V, arma::mat theta, int Lnm)//output of form_Fv(*)
{
	arma::mat L = arma::zeros(Lnm - 1, Lnm - 1);
	const double pi = 3.14159265358979323846;
	for (int i = 1; i < Lnm; ++i)
	{
		long double R = 0;
		for (int j = 1; j < Lnm; ++j)
		{
			if (i != j)
			{
				L(i - 1, j - 1) = V(i, 0)*(real(Y(i, j))*sin((theta(0, i) - theta(0, j))*pi / 180) - imag(Y(i, j))*cos((theta(0, i) - theta(0, j))*pi / 180));
			}
			else
			{
				for (int m = 0; m < Lnm; ++m)
				{
					if (m != i)
					{
						R = R + V(m, 0)*(real(Y(i, m))*sin((theta(0, i) - theta(0, m))*pi / 180) - imag(Y(i, m))*cos((theta(0, i) - theta(0, m))*pi / 180));
					}
				}
				L(i - 1, j - 1) = -2 * V(i, 0)*imag(Y(i, j)) + R;
			}
		}
	}
	return L;
}