//# define ARMA_DONT_USE_WRAPPER
# include <armadillo>
# include <iostream>
# include <math.h>
# include "fun_return.h"
#include "load_system_data.h"

using namespace std;
using namespace arma;

vrmove VR_move(mat Vpolar, mat Dl, mat v_lim, int Lvp)
{
	vrmove vr_move;
	double lb, ub, epsi;
	epsi = 0.0001;
	ub = v_lim(0) + epsi;
	lb = v_lim(1) - epsi;
	cout << ub << endl << lb << endl;

	Mat<int> VR1_zone, VR2_zone, VR3_zone, VR4_zone;
	std::ifstream file1("matlab_data_VR1.txt");
	VR1_zone.load(file1);
	std::ifstream file2("matlab_data_VR2.txt");
	VR2_zone.load(file2);
	std::ifstream file3("matlab_data_VR3.txt");
	VR3_zone.load(file3);
	std::ifstream file4("matlab_data_VR4.txt");
	VR4_zone.load(file4);

	//cout << "VR1_zone: \n" << VR1_zone << endl;
	//mat VR1_edge, VR2_edge, VR3_edge, VR4_edge;
	//VR1_edge << 0 << 1 << endr;
	//VR2_edge << 12 << 13 << endr;
	//VR3_edge << 30 << 31 << endr;
	//VR4_edge << 90 << 118 << endr;

	mat VR1, VR2, VR3, VR4;
	VR1 << 0 << 0 << 0 << endr;
	VR2 << 0 << 0 << 0 << endr;
	VR3 << 0 << 0 << 0 << endr;
	VR4 << 0 << 0 << 0 << endr;

	int j = 0;
	for (int i = 0; i < Lvp; ++i)
	{
		if (Vpolar(i, 0) == 0)
		{
			Vpolar(i, 0) = NAN;
		}
		if (Vpolar(i, 2) == 0)
		{
			Vpolar(i, 2) = NAN;
		}
		if (Vpolar(i, 4) == 0)
		{
			Vpolar(i, 4) = NAN;
		}
	}

	int phase;
	//mat Location_Vh;
	for (phase = 0; phase < 3; phase++)
	{
		//int j = 0;
		//for (int i = 0; i < Lvp; ++i)
		//{
		//	if (Vpolar(i, 2 * phase) > ub)
		//	{
		//		Location_Vh(j, phase) = i;
		//		j++;
		//	}
		//}

		uvec Location_Vh = find(Vpolar.col(2 * phase) > ub);
		//cout << "Location_Vh" << Location_Vh << endl;

		int down_VR1, down_VR2, down_VR3, down_VR4;
		down_VR1 = 0;
		down_VR2 = 0;
		down_VR3 = 0;
		down_VR4 = 0;

		if (!Location_Vh.is_empty())
		{
			for (int k = 0; k < Location_Vh.n_rows; ++k)
			{
				uvec temph1 = find(VR1_zone == Location_Vh(k));
				uvec temph2 = find(VR2_zone == Location_Vh(k));
				uvec temph3 = find(VR3_zone == Location_Vh(k));
				uvec temph4 = find(VR4_zone == Location_Vh(k));
				if (!temph1.is_empty())
				{
					down_VR1 = 1;

					if (!temph2.is_empty())
					{
						down_VR2 = 1;
					}
					else if (!temph3.is_empty())
					{
						down_VR3 = 1;
					}
					else if (!temph4.is_empty())
					{
						down_VR4 = 1;
					}
					else
					{
					}
				}
				else
				{
					cout << "The node that has voltage violdation doesnot belong to any VR zones: " << Location_Vh(k) << endl;
				}
			}
		}

		uvec Location_Vl = find(Vpolar.col(2 * phase) < lb);
		//cout << "Location_Vl" << Location_Vl << endl;

		int up_VR1, up_VR2, up_VR3, up_VR4;
		up_VR1 = 0;
		up_VR2 = 0;
		up_VR3 = 0;
		up_VR4 = 0;

		if (!Location_Vl.is_empty())
		{
			for (int k = 0; k < Location_Vl.n_rows; ++k)
			{
				uvec templ1 = find(VR1_zone == Location_Vl(k));
				uvec templ2 = find(VR2_zone == Location_Vl(k));
				uvec templ3 = find(VR3_zone == Location_Vl(k));
				uvec templ4 = find(VR4_zone == Location_Vl(k));
				if (!templ1.is_empty())
				{
					up_VR1 = 1;

					if (!templ2.is_empty())
					{
						up_VR2 = 1;
					}
					else if (!templ3.is_empty())
					{
						up_VR3 = 1;
					}
					else if (!templ4.is_empty())
					{
						up_VR4 = 1;
					}
					else
					{
					}
					//cout << "found location!" << endl;
				}
				else
				{
					cout << "The node that has voltage violdation doesnot belong to any VR zones: " << Location_Vl(k) << endl;
				}
			}
		}

		if (up_VR1 == 1 && down_VR1 == 1)
		{
			cout << "Large Voltage Variation in Zone of VR1 on Phase: " << phase << endl;
		}
		else
		{
			VR1(phase) = up_VR1 - down_VR1;
		}

		if (up_VR2 == 1 && down_VR2 == 1)
		{
			cout << "Large Voltage Variation in Zone of VR2 on Phase: " << phase << endl;
		}
		else
		{
			VR2(phase) = up_VR2 - down_VR2;
		}

		if (up_VR3 == 1 && down_VR3 == 1)
		{
			cout << "Large Voltage Variation in Zone of VR3 on Phase: " << phase << endl;
		}
		else
		{
			VR3(phase) = up_VR3 - down_VR3;
		}

		if (up_VR4 == 1 && down_VR4 == 1)
		{
			cout << "Large Voltage Variation in Zone of VR4 on Phase: " << phase << endl;
		}
		else
		{
			VR4(phase) = up_VR4 - down_VR4;
		}
	}
	cout << "Exiting VR_move.cpp" << endl;
	vr_move.VR1 = VR1;
	vr_move.VR2 = VR2;
	vr_move.VR3 = VR3;
	vr_move.VR4 = VR4;
	return vr_move;
}