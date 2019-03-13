#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <string>
#include <fstream>
#include <boost/array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/numeric/odeint.hpp>

#define NUMPOINTS 1500

typedef boost::array<double, 9> cell_state;
typedef boost::unordered_map<int, boost::array<double, 9>> population;

void linspace(boost::array<double, NUMPOINTS> &SNAIL, double start_val, double end_val)
{
	if(start_val <= end_val)
	{
		for(int i = 0; i < NUMPOINTS; i++)
		{
			SNAIL[i] = start_val + i*(end_val - start_val) / (double) NUMPOINTS;
		}
	}
	else
	{
		for(int i = 0; i < NUMPOINTS; i++)
		{
			SNAIL[i] = start_val - i*(start_val - end_val) / (double) NUMPOINTS;
		}
	}
}

double nchoosek(int n, int k)
{
	double result = 1.0;
	double result0 = 1.0;

	for(int i = 0; i < k; i++)
	{
		result *= (n - i);
	}
	for(int i = 1; i <= k; i++)
	{
		result0 *= i;
	}

	return(result / result0);
}

void SNAIL_ZEB_miR200_miR34_GRHL2_system(const cell_state &x, cell_state &dxdt, double t)
{
	double g_miR34 = 1.35e3;
	double g_mSNAIL = 90.0;
	double g_SNAIL = 0.1e3;
	double g_miR200 = 2.1e3;
	double g_mZEB = 11.0;
	double g_ZEB = 0.1e3;
	double g_mG = 22.0;
	double g_G = 200.0;

	double k_miR34 = 0.05;
	double k_mSNAIL = 0.5;
	double k_SNAIL = 0.125;
	double k_miR200 = 0.05;
	double k_mZEB = 0.5;
	double k_ZEB = 0.1;
	double k_mG = 0.5;
	double k_G = 0.1;

	double t_miR34_SNAIL = 300e3;
	double t_mSNAIL_SNAIL = 200e3;
	double t_miR34_ZEB = 600e3;
	double t_miR34 = 10e3;
	double t_mSNAIL_I = 50e3;
	double t_miR200_ZEB = 220e3;
	double t_miR200_SNAIL = 180e3;
	double t_mZEB_ZEB = 25e3;
	double t_mZEB_SNAIL = 180e3;
	double t_miR200 = 10e3;
	double t_mZEB_G = 25e3;
	double t_mG_ZEB = 10e3;

	int n_miR34_SNAIL = 1;
	int n_miR34_ZEB = 1;
	int n_mSNAIL_SNAIL = 1;
	int n_mSNAIL_I = 1;
	int n_miR200_ZEB = 3;
	int n_miR200_SNAIL = 2;
	int n_mZEB_ZEB = 2;
	int n_mZEB_SNAIL = 2;
	int n_mZEB_G = 1;
	int n_mG_ZEB = 3;

	double l_miR34_SNAIL = 0.1;
	double l_mSNAIL_SNAIL = 0.1;
	double l_miR34_ZEB = 0.2;
	double l_mSNAIL_I = 10.0;
	double l_miR200_ZEB = 0.1;
	double l_miR200_SNAIL = 0.1;
	double l_mZEB_ZEB = 7.5;
	double l_mZEB_SNAIL = 10.0;
	double l_mZEB_G = 0.65;
	double l_mG_ZEB = 0.5;

	boost::array<double, 7> L = {1.0, 0.6, 0.3, 0.1, 0.05, 0.05, 0.05};
	boost::array<double, 7> gamma_mRNA = {0.0, 0.04, 0.2, 1.0, 1.0, 1.0, 1.0};
	boost::array<double, 7> gamma_miRNA = {0.0, 0.005, 0.05, 0.5, 0.5, 0.5, 0.5};

	double degrad_miR200 = 0.0;
	double degrad_mZEB = 0.0;
	double trans_mZEB = 0.0;

	double fac_miR200 = 0.0;

	for(int i = 0; i < 7; i++)
	{
		fac_miR200 = (std::pow(x[0] / t_miR200, i)) / (std::pow(1 + (x[0] / t_miR200), 6));
		degrad_miR200 += gamma_miRNA[i]*nchoosek(6, i)*i*fac_miR200;
		degrad_mZEB += gamma_mRNA[i]*nchoosek(6, i)*fac_miR200;
		trans_mZEB += L[i]*nchoosek(6, i)*fac_miR200;
	}

	double degrad_miR34 = 0.0;
	double degrad_mSNAIL = 0.0;
	double trans_mSNAIL = 0.0;

	double fac_miR34 = 0.0;

	for(int i = 0; i < 3; i++)
	{
		fac_miR34 = (std::pow(x[5] / t_miR34, i)) / (std::pow(1 + (x[5] / t_miR34), 2));
		degrad_miR34 += gamma_miRNA[i]*nchoosek(2, i)*i*fac_miR34;
		degrad_mSNAIL += gamma_mRNA[i]*nchoosek(2, i)*fac_miR34;
		trans_mSNAIL += L[i]*nchoosek(2, i)*fac_miR34;
	}

	double H_miR200_ZEB = 1 / (1 + std::pow(x[2] / t_miR200_ZEB, n_miR200_ZEB));
	H_miR200_ZEB = H_miR200_ZEB + l_miR200_ZEB*(1 - H_miR200_ZEB);

	double H_miR200_SNAIL = 1 / (1 + std::pow(x[3] / t_miR200_SNAIL, n_miR200_SNAIL));
	H_miR200_SNAIL = H_miR200_SNAIL + l_miR200_SNAIL*(1 - H_miR200_SNAIL);

	double H_mZEB_ZEB = 1 / (1 + std::pow(x[2] / t_mZEB_ZEB, n_mZEB_ZEB));
	H_mZEB_ZEB = H_mZEB_ZEB + l_mZEB_ZEB*(1 - H_mZEB_ZEB);

	double H_mZEB_SNAIL = 1 / (1 + std::pow(x[3] / t_mZEB_SNAIL, n_mZEB_SNAIL));
	H_mZEB_SNAIL = H_mZEB_SNAIL + l_mZEB_SNAIL*(1 - H_mZEB_SNAIL);

	double H_miR34_SNAIL = 1 / (1 + std::pow(x[3] / t_miR34_SNAIL, n_miR34_SNAIL));
	H_miR34_SNAIL = H_miR34_SNAIL + l_miR34_SNAIL*(1 - H_miR34_SNAIL);

	double H_miR34_ZEB = 1 / (1 + std::pow(x[2] / t_miR34_ZEB, n_miR34_ZEB));
	H_miR34_ZEB = H_miR34_ZEB + l_miR34_ZEB*(1 - H_miR34_ZEB);

	double H_mSNAIL_SNAIL = 1 / (1 + std::pow(x[3] / t_mSNAIL_SNAIL, n_mSNAIL_SNAIL));
	H_mSNAIL_SNAIL = H_mSNAIL_SNAIL + l_mSNAIL_SNAIL*(1 - H_mSNAIL_SNAIL);

	double H_mSNAIL_I = 1 / (1 + std::pow(x[6] / t_mSNAIL_I, n_mSNAIL_I));
	H_mSNAIL_I = H_mSNAIL_I + l_mSNAIL_I*(1 - H_mSNAIL_I);

	double H_mZEB_G = 1 / (1 + std::pow(x[8] / t_mZEB_G, n_mZEB_G));
	H_mZEB_G = H_mZEB_G + l_mZEB_G*(1 - H_mZEB_G);

	double H_mG_ZEB = 1 / (1 + std::pow(x[2] / t_mG_ZEB, n_mG_ZEB));
	H_mG_ZEB = H_mG_ZEB + l_mG_ZEB*(1 - H_mG_ZEB);

	dxdt[0] = g_miR200*H_miR200_ZEB*H_miR200_SNAIL - x[1]*degrad_miR200 - k_miR200*x[0];
	dxdt[1] = g_mZEB*H_mZEB_ZEB*H_mZEB_SNAIL*H_mZEB_G - x[1]*degrad_mZEB - k_mZEB*x[1];
	dxdt[2] = g_ZEB*x[1]*trans_mZEB - k_ZEB*x[2];
	dxdt[3] = g_SNAIL*x[4]*trans_mSNAIL - k_SNAIL*x[3];
	dxdt[4] = g_mSNAIL*H_mSNAIL_I*H_mSNAIL_SNAIL - x[4]*degrad_mSNAIL - k_mSNAIL*x[4];
	dxdt[5] = g_miR34*H_miR34_ZEB*H_miR34_SNAIL - x[4]*degrad_miR34 - k_miR34*x[5];
	dxdt[6] = 0.0;
	dxdt[7] = g_mG*H_mG_ZEB - k_mG*x[7];
	dxdt[8] = g_G*x[7] - k_G*x[8];
}

void get_bifurcation(double start_val, double end_val)
{
	boost::array <double, NUMPOINTS> I;
	boost::array <double, NUMPOINTS> mZEB;

	linspace(I, start_val, end_val);
	cell_state single_cell;
	single_cell = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	for(int i = 0; i < NUMPOINTS; i++)
	{
		single_cell[6] = I[i];
		boost::numeric::odeint::integrate(SNAIL_ZEB_miR200_miR34_GRHL2_system, single_cell, 0.0, 7500.0, 0.1);
		mZEB[i] = single_cell[1];
	}
	for(int i = 0; i < NUMPOINTS; i++)
	{
		std::cout << I[i] << " " << mZEB[i] << "\n";
	}
}

int main(int argc, char *argv[])
{
	get_bifurcation(20e3, 220e3);
	get_bifurcation(220e3, 20e3);
	get_bifurcation(130e3, 20e3);
	return(0);
}