#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <string>
#include <fstream>
#include <boost/array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/numeric/odeint.hpp>
//#include <mpi.h>

#define NUMNODES 7

typedef boost::array<double, NUMNODES> cell_state;
typedef boost::unordered_map<int, boost::array<double, NUMNODES>> population;

std::mt19937 generator;

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

void SNAIL_ZEB_miR200_miR34_system(const cell_state &x, cell_state &dxdt, double t)
{
	double g_miR34 = 1.35e3;
	double g_mSNAIL = 90.0;
	double g_SNAIL = 0.1e3;
	double g_miR200 = 2.1e3;
	double g_mZEB = 11.0;
	double g_ZEB = 0.1e3;

	double k_miR34 = 0.05;
	double k_mSNAIL = 0.5;
	double k_SNAIL = 0.125;
	double k_miR200 = 0.05;
	double k_mZEB = 0.5;
	double k_ZEB = 0.1;

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

	int n_miR34_SNAIL = 1;
	int n_miR34_ZEB = 1;
	int n_mSNAIL_SNAIL = 1;
	int n_mSNAIL_I = 1;
	int n_miR200_ZEB = 3;
	int n_miR200_SNAIL = 2;
	int n_mZEB_ZEB = 2;
	int n_mZEB_SNAIL = 2;

	double l_miR34_SNAIL = 0.1;
	double l_mSNAIL_SNAIL = 0.1;
	double l_miR34_ZEB = 0.2;
	double l_mSNAIL_I = 10.0;
	double l_miR200_ZEB = 0.1;
	double l_miR200_SNAIL = 0.1;
	double l_mZEB_ZEB = 7.5;
	double l_mZEB_SNAIL = 10.0;

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

	dxdt[0] = g_miR200*H_miR200_ZEB*H_miR200_SNAIL - x[1]*degrad_miR200 - k_miR200*x[0];
	dxdt[1] = g_mZEB*H_mZEB_ZEB*H_mZEB_SNAIL - x[1]*degrad_mZEB - k_mZEB*x[1];
	dxdt[2] = g_ZEB*x[1]*trans_mZEB - k_ZEB*x[2];
	dxdt[3] = g_SNAIL*x[4]*trans_mSNAIL - k_SNAIL*x[3];
	dxdt[4] = g_mSNAIL*H_mSNAIL_I*H_mSNAIL_SNAIL - x[4]*degrad_mSNAIL - k_mSNAIL*x[4];
	dxdt[5] = g_miR34*H_miR34_ZEB*H_miR34_SNAIL - x[4]*degrad_miR34 - k_miR34*x[5];
	dxdt[6] = 0.0;
}

void initialize_Signal_lognormal(population &P, int N)
{
	std::normal_distribution <> dist{0.0, 1.0};

	double CV2 = 1.0;
	double SD = std::sqrt(std::log(CV2*CV2 + 1.0));
	double M = 20e3;

	for(int i = 0; i < N; i++)
	{
		P[i] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		P[i][6] = std::exp(std::log(M) + dist(generator)*SD);
	}

	for(int i = 0; i < N; i++)
	{
		boost::numeric::odeint::integrate(SNAIL_ZEB_miR200_miR34_system, P[i], 0.0, 7500.0, 0.1);
	}
}

void get_phenotypes(population &P, boost::unordered_map <int, int> &phenotype)
{
	double x1 = 30.1, x2 = 66.5, y1 = 719.537, y2 = 477.792;
	double a1 = y2 - y1, b1 = x1 - x2, c1 = x2*y1 - y2*x1;

	double u1 = 45.65, u2 = 59.9, v1 = 252.173, v2 = 89.4267;
	double a2 = v2 - v1, b2 = u1 - u2, c2 = u2*v1 - v2*u1;

	double x, y, fac1, fac2;
	int state = -1;

	for(int i = 0; i < P.size(); i++)
	{
		x = P[i][6] / 1e3;
		y = P[i][1];
		state = -1;

		if(x <= x1)
		{
			state = 0;
		}
		else if(x > x2)
		{
			state = 2;
		}
		else
		{
			fac1 = (a1*x + b1*y + c1) / b1;
			fac2 = (a2*x + b2*y + c2) / b2;
			if(fac1 >= 0 && x > x1 && x < x2)
			{
				state = 2;
			}
			else if(fac1 < 0 && x < u1)
			{
				state = 0;
			}
			else if(fac2 < 0 && x < u2)
			{
				state = 0;
			}
			else if(fac2 >= 0 && fac1 < 0 && x > u1 && x < u2)
			{
				state = 1;
			}
			else if(fac1 < 0 && x > u2 && x < x2)
			{
				state = 1;
			}
		}
		if(state == -1)
		{
			std::cout << "Error in phenotype assignment." << "\n";
		}
		else
		{
			phenotype[i] = state;
		}
	}
}

int get_sample(boost::unordered_map <int, int> &phenotype, int pid)
{
	std::vector<int> V;
	int count = 0;
	for(int i = 0; i < phenotype.size(); i++)
	{
		if(phenotype[i] == pid)
		{
			V.push_back(i);
			count += 1;
		}
	}

	std::uniform_int_distribution<int> distribution(0, count - 1);

	return(V[distribution(generator)]);
}

void add_asym_noise_normal(population &P, int rep_id, int new_id, int rep_phenotype, double eta)
{
	if(rep_phenotype == 1)
	{
		P[new_id][5] = 0.0;
	}

	std::normal_distribution <> distribution{0, 1};

        P[rep_id][6] = P[rep_id][6] + distribution(generator)*eta;
        P[new_id][6] = P[new_id][6] + distribution(generator)*eta;

        if(P[rep_id][6] < 0.0)
        {
                P[rep_id][6] = 0.0;
        }
        if(P[new_id][6] < 0.0)
        {
                P[new_id][6] = 0.0;
        }
}

void simulate_normal(population &P, double end_time, double eta, int eta_id, int rank, boost::array <double, 3> &GR)
{
	boost::unordered_map <int, double> update_time;
	boost::unordered_map <int, int> phenotype;
	boost::array<double, 3> count;
	boost::unordered_map <int, boost::array<int, 7>> TrackDiv;
	int DivIndex = 0;

	for(int i = 0; i < P.size(); i++)
	{
		update_time[i] = 0.0;
	}

	int last_index = P.size();

	get_phenotypes(P, phenotype);

	boost::array <double, 3> growth;
	growth[0] = std::log(2) / GR[0];
	growth[1] = std::log(2) / GR[1];
	growth[2] = std::log(2) / GR[2];
	boost::array<double, 3> r0 = {growth[0], growth[1], growth[2]};
	boost::array<double, 3> d0 = {growth[0] / 10.0, growth[1] / 10.0, growth[2] / 10.0};
	boost::array<double, 6> rates = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	int K = 10000;

	for(int i = 0; i < P.size(); i++)
	{
		if(phenotype[i] == 0)
		{
			if(P.size() < K)
			{
				rates[0] += r0[0]*(1 - float(P.size()) / K);
			}
			rates[3] += d0[0];
		}
		else if(phenotype[i] == 1)
		{
			if(P.size() < K)
			{
				rates[1] += r0[1]*(1 - float(P.size()) / K);
			}
			rates[4] += d0[1];
		}
		else
		{
			if(P.size() < K)
			{
				rates[2] += r0[2]*(1 - float(P.size()) / K);
			}
			rates[5] += d0[2];
		}
	}

	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	double t = 0.0, a0, dt, p0, p1;
	double sum_prev, sum_new;
	double last_updated = 0.0;
	int event_id, rep_id, death_id;
	std::vector<int> ToErase;

	count = {0.0, 0.0, 0.0};
	for(int i = 0; i < P.size(); i++)
	{
		count[phenotype[i]] += 1.0;
	}
	for(int i = 0; i < 3; i++)
	{
		count[i] = double(count[i]) / P.size();
	}

	std::ofstream count_file;
	std::ofstream TrackFile;

	count_file.open("outputfiles/RUN" + std::to_string(eta_id) + "/count_file_" + std::to_string(eta_id) + '_' + std::to_string(rank) + ".log");
	count_file << t / 24.0 << " " << count[0] << " " << count[1] << " " << count[2] << " " << P.size() << "\n";
	TrackFile.open("outputfiles/RUN" + std::to_string(eta_id) + "/TrackFile_" + std::to_string(eta_id) + '_' + std::to_string(rank) + ".log");

	while(t <= end_time)
	{
		a0 = 0.0;

		for(int i = 0; i < 6; i++)
		{
			a0 += rates[i];
		}

		p0 = distribution(generator);
		dt = (1.0 / a0)*std::log(1 / p0);
		t += dt;

		p1 = distribution(generator);
		sum_prev = 0.0;
		sum_new = 0.0;
		event_id = -1;

		for(int i = 0; i < 6; i++)
		{
			sum_new = sum_prev + rates[i] / a0;
			if(p1 >= sum_prev && p1 < sum_new)
			{
				event_id = i;
				break;
			}
			sum_prev = sum_new;
		}

		if(event_id == 0)
		{
			rep_id = get_sample(phenotype, 0);
		}
		else if(event_id == 1)
		{
			rep_id = get_sample(phenotype, 1);
		}
		else if(event_id == 2)
		{
			rep_id = get_sample(phenotype, 2);
		}
		else if(event_id == 3)
		{
			death_id = get_sample(phenotype, 0);
		}
		else if(event_id == 4)
		{
			death_id = get_sample(phenotype, 1);
		}
		else if(event_id == 5)
		{
			death_id = get_sample(phenotype, 2);
		}

		if(event_id < 3)
		{
			boost::numeric::odeint::integrate(SNAIL_ZEB_miR200_miR34_system, P[rep_id], 0.0, t - update_time[rep_id], 0.1);
			get_phenotypes(P, phenotype);
			for(auto p: TrackDiv)
			{
				for(int k = 1; k < 3; k++)
				{
					if(p.second[k] == rep_id && p.second[k + 4] == 0)
					{
						TrackDiv[p.first][k + 2] = phenotype[rep_id];
						TrackDiv[p.first][k + 4] = 1;
					}
				}
			}
			P[last_index] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			for(int j = 0; j < NUMNODES; j++)
			{
				P[last_index][j] = P[rep_id][j];
			}
			phenotype[last_index] = phenotype[rep_id];
			add_asym_noise_normal(P, rep_id, last_index, phenotype[rep_id], eta);
			update_time[rep_id] = t;
			update_time[last_index] = t;
			last_index += 1;
			TrackDiv[DivIndex] = {phenotype[rep_id], rep_id, last_index - 1, phenotype[rep_id], phenotype[last_index - 1], 0, 0};
			DivIndex += 1;
		}
		else
		{
			boost::numeric::odeint::integrate(SNAIL_ZEB_miR200_miR34_system, P[death_id], 0.0, t - update_time[death_id], 0.1);
			get_phenotypes(P, phenotype);
			for(auto p: TrackDiv)
			{
				for(int k = 1; k < 3; k++)
				{
					if(p.second[k] == death_id && p.second[k + 4] == 0)
					{
						TrackDiv[p.first][k + 2] = phenotype[death_id];
						TrackDiv[p.first][k + 4] = 1;
					}
					if(p.second[k] == last_index - 1 && p.second[k + 4] == 0)
					{
						TrackDiv[p.first][k] = death_id;
					}
				}
			}
			for(int j = 0; j < NUMNODES; j++)
			{
				P[death_id][j] = P[last_index - 1][j];
			}
			update_time[death_id] = update_time[last_index - 1];
			phenotype[death_id] = phenotype[last_index - 1];
			P.erase(last_index - 1);
			update_time.erase(last_index - 1);
			phenotype.erase(last_index - 1);
			last_index -= 1;
		}

		if(true)
		{
			for(int i = 0; i < P.size(); i++)
			{
				if(t - update_time[i] > 0.0)
				{
					boost::numeric::odeint::integrate(SNAIL_ZEB_miR200_miR34_system, P[i], 0.0, t - update_time[i], 0.1);
				}
				update_time[i] = t;
			}
		}

		get_phenotypes(P, phenotype);

		rates = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

		for(int i = 0; i < P.size(); i++)
		{
			if(phenotype[i] == 0)
			{
				if(P.size() < K)
				{
					rates[0] += r0[0]*(1 - float(P.size()) / K);
				}
				rates[3] += d0[0];
			}
			else if(phenotype[i] == 1)
			{
				if(P.size() < K)
				{
					rates[1] += r0[1]*(1 - float(P.size()) / K);
				}
				rates[4] += d0[1];
			}
			else
			{
				if(P.size() < K)
				{
					rates[2] += r0[2]*(1 - float(P.size()) / K);
				}
				rates[5] += d0[2];
			}
		}

		if(true)
		{
			count = {0.0, 0.0, 0.0};
			for(int i = 0; i < P.size(); i++)
			{
				count[phenotype[i]] += 1.0;
			}
			for(int i = 0; i < 3; i++)
			{
				count[i] = double(count[i]) / P.size();
			}
			count_file << t / 24.0 << " " << count[0] << " " << count[1] << " " << count[2] << " " << P.size() << "\n";
			last_updated = t;
		}
		for(auto p: TrackDiv)
		{
			if(p.second[5] == 1 && p.second[6] == 1)
			{
				TrackFile << p.second[0] << " " << p.second[3] << " " << p.second[4] << "\n";
				ToErase.push_back(p.first);
			}
		}
		for(int j = 0; j < ToErase.size(); j++)
		{
			TrackDiv.erase(ToErase[j]);
		}
	}

	TrackFile.close();
	count_file.close();
}

int main(int argc, char *argv[])
{
	//MPI_Init(NULL, NULL);
	int world_size;
	//MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank = 0;
	//MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	boost::array <double, 3> GR;

	double eta = std::stod(argv[1])*1e3;
	int eta_id = std::stoi(argv[2]);
	GR[0] = std::stod(argv[3]);
	GR[1] = std::stod(argv[4]);
	GR[2] = std::stod(argv[5]);

	generator = std::mt19937(std::time(NULL) + world_rank + eta_id);

	population P;
	int N = 500;

	initialize_Signal_lognormal(P, N);
	simulate_normal(P, 24.0*56, eta, eta_id, world_rank, GR);

	//MPI_Finalize();

	return(0);
}
