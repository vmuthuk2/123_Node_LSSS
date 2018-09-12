#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>

int main()
{
	arma::Mat<double> Dl;
	arma::Mat<double> R;
	arma::Mat<double> X;
	arma::cx_mat<double> Z;

	std::ifstream file1("matlab_data_Dl.txt");
	Dl.load(file);
	std::ifstream file1("matlab_data_R.txt");
	R.load(file);
	std::ifstream file1("matlab_data_X.txt");
	X.load(file);
	Z = arma::cx_mat(R, X);
	return 0;
}