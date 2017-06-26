#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<algorithm>
#include<vector>
#include<list>
#include<queue>
#include<stack>
#include<fstream>
#include<iomanip>
#include<ctime>
#include<cmath>
#include<direct.h>
//#include"mersenne.h"
#include "MT.h"
#include"wfg.h"
#include"moead.h"
#include <chrono>

using namespace std;

#define X_MIN 0.0
#define X_MAX 1.0
#define EPS 0.00000000000001



void SBXcrossover(double *par_x1, double *par_x2, double **off_x, int off_p_n, int item, double CO, double SBX_ETA){
	int off_p_n1 = off_p_n;
	int off_p_n2 = off_p_n + 1;

	double x1 = 0.0;
	double x2 = 0.0;
	double beta1 = 1.0;
	double beta2 = 1.0;
	double alpha1 = 2.0;
	double alpha2 = 2.0;
	double beta_q_1 = 1.0;
	double beta_q_2 = 1.0;
	double c1, c2;
	double u = nextDoubleIE();

	if (nextDoubleIE() < CO){
		for (int n = 0; n < item; n++){
			if (nextDoubleIE() < 0.5){
				if (abs(par_x2[n] - par_x1[n]) > EPS){
					x1 = 0.0;
					x2 = 0.0;
					beta1 = 1.0;
					beta2 = 1.0;
					alpha1 = 2.0;
					alpha2 = 2.0;
					beta_q_1 = 1.0;
					beta_q_2 = 1.0;

					u = nextDoubleIE();
					if (par_x1[n] > par_x2[n]){
						x2 = par_x1[n];
						x1 = par_x2[n];
					}
					else if (par_x1[n] < par_x2[n]){
						x2 = par_x2[n];
						x1 = par_x1[n];
					}


					beta1 += 2.0 * (x1 - X_MIN) / (x2 - x1);
					beta2 += 2.0 * (X_MAX - x2) / (x2 - x1);

					alpha1 -= pow(beta1, -1.0 * (SBX_ETA + 1.0));
					alpha2 -= pow(beta2, -1.0 * (SBX_ETA + 1.0));

					if (u <= 1.0 / alpha1){
						beta_q_1 = pow(u * alpha1, (1.0 / (SBX_ETA + 1.0)));
					}
					else{
						beta_q_1 = pow(1.0 / (2.0 - u * alpha1), (1.0 / (SBX_ETA + 1.0)));
					}

					if (u <= 1.0 / alpha2){
						beta_q_2 = pow(u * alpha2, (1.0 / (SBX_ETA + 1.0)));
					}
					else{
						beta_q_2 = pow(1.0 / (2.0 - u * alpha2), (1.0 / (SBX_ETA + 1.0)));
					}

					c1 = 0.5 *((x1 + x2) - beta_q_1 * (x2 - x1));
					c2 = 0.5 *((x1 + x2) + beta_q_2 * (x2 - x1));

					if (c1 < X_MIN) c1 = X_MIN;
					if (c1 > X_MAX) c1 = X_MAX;

					if (c2 < X_MIN) c2 = X_MIN;
					if (c2 > X_MAX) c2 = X_MAX;
					
					
					//improvement
					/*off_x[off_p_n1][n] = c1;
					off_x[off_p_n2][n] = c2;
					*/
					
				
					if (nextDoubleIE() < 0.5){
						off_x[off_p_n1][n] = c1;
						off_x[off_p_n2][n] = c2;

					}
					else{
						off_x[off_p_n1][n] = c2;
						off_x[off_p_n2][n] = c1;
					}
					//improvement

					
				}
				else{
					off_x[off_p_n1][n] = par_x1[n];
					off_x[off_p_n2][n] = par_x2[n];
				}
			}
			else{
				off_x[off_p_n1][n] = par_x1[n];
				off_x[off_p_n2][n] = par_x2[n];
			}
		}

	}
	else{
		for (int n = 0; n < item; n++){
			off_x[off_p_n1][n] = par_x1[n];
			off_x[off_p_n2][n] = par_x2[n];
		}

	}


}

void polynomial_mutation(double *off_x, int item, double MUT, double MUT_ETA){

	double r = 0.0;
	double shita1 = 0.0;
	double shita2 = 0.0;
	double shitaq = 0.0;
	for (int i = 0; i < item; i++){
		r = nextDoubleIE();
		if (r <= MUT){
			shita1 = (off_x[i] - X_MIN) / (X_MAX - X_MIN);
			shita2 = (X_MAX - off_x[i]) / (X_MAX - X_MIN);

			r = nextDoubleIE();

			if (r <= 0.5){
				shitaq = pow(((2.0 * r) + (1.0 - 2.0 * r) * pow(1.0 - shita1, MUT_ETA + 1.0)), 1.0 / (MUT_ETA + 1.0)) - 1.0;

				off_x[i] += shitaq * (X_MAX - X_MIN);

			}
			else{
				shitaq = 1.0 - pow(2.0 * (1.0 - r) + 2.0 * (r - 0.5) * pow(1.0 - shita2, MUT_ETA + 1.0), 1.0 / (MUT_ETA + 1.0));
				off_x[i] += shitaq * (X_MAX - X_MIN);
			}
			if (off_x[i] < X_MIN){
				off_x[i] = X_MIN;
			}
			if (off_x[i] > X_MAX){
				off_x[i] = X_MAX;
			}
		}

	}
}


//continuous optimization
void make_new_pop(double *x1, double *x2, double off_x[], int item, double MUT, double CO,double MUT_ETA, double SBX_ETA){

	double **temp_off_x = new double*[2];

	temp_off_x[0] = new double[item];
	temp_off_x[1] = new double[item];

	SBXcrossover(x1, x2, temp_off_x, 0, item, CO, SBX_ETA);
	if (nextDoubleIE() < 0.5){
		for (int i = 0; i < item; i++){
			off_x[i] = temp_off_x[0][i];
		}
	}
	else{
		for (int i = 0; i < item; i++){
			off_x[i] = temp_off_x[1][i];
		}
	}

	polynomial_mutation(off_x, item, MUT, MUT_ETA);

	for (int i = 0; i < 2; i++){
		delete[] temp_off_x[i];
	}
	delete[] temp_off_x;

}


void make_new_pop(int *x1, int *x2, int off_x[], int item, double MUT, double CO){

	int *temp = new int[item];
	int **x = new int *[2];
	x[0] = new int[item];
	x[1] = new int[item];

	for (int i = 0; i < item; i++){
		x[0][i] = x1[i];
	}
	for (int i = 0; i < item; i++){
		x[1][i] = x2[i];
	}

	if (nextDoubleIE() < CO){
		for (int i = 0; i < item; i++){
			temp[i] = genrand_int31() % 2;
		}
		/*
		for (int i = 0; i < item; i++){
		off_x[i] = x[temp[i]][i];
		if (nextDoubleIE() <= MUT){
		off_x[i] = (off_x[i] + 1) % 2;
		}
		}*/
	}
	else{
		int k = genrand_int31() % 2;
		for (int i = 0; i < item; i++){
			temp[i] = k;
			/*	off_x[i] = x[temp[i]][i]; */
		}
	}

	for (int i = 0; i < item; i++){
		off_x[i] = x[temp[i]][i];
		if (nextDoubleIE() < MUT){
			off_x[i] = (off_x[i] + 1) % 2;
		}
	}
	delete[] temp;
	for (int n = 0; n < 2; n++){
		delete[] x[n];
	}
	delete[] x;
}

void make_new_pop_20160301(int *x1, int *x2, int off_x[], int item, double MUT, double CO, int ob){

	int *temp = new int[item];
	int **x = new int *[2];
	x[0] = new int[item];
	x[1] = new int[item];

	for (int i = 0; i < item; i++){
		x[0][i] = x1[i];
	}
	for (int i = 0; i < item; i++){
		x[1][i] = x2[i];
	}

	if (nextDoubleIE() < CO){
		for (int i = 0; i < item; i++){
			temp[i] = genrand_int31() % 2;
		}
		/*
		for (int i = 0; i < item; i++){
		off_x[i] = x[temp[i]][i];
		if (nextDoubleIE() <= MUT){
		off_x[i] = (off_x[i] + 1) % 2;
		}
		}*/
	}
	else{
		int k = genrand_int31() % 2;
		for (int i = 0; i < item; i++){
			temp[i] = k;
			/*	off_x[i] = x[temp[i]][i]; */
		}
	}

	for (int i = 0; i < item; i++){
		off_x[i] = x[temp[i]][i];
		if (nextDoubleIE() < MUT){
			int t = off_x[i];
			do{
				off_x[i] = genrand_int31() % ob + 1;
			} while (t == off_x[i]);
		}
	}
	delete[] temp;
	for (int n = 0; n < 2; n++){
		delete[] x[n];
	}
	delete[] x;
}


bool dominate(double *fita, double *fitb, int ob, int SIGN){

	for (int o = 0; o < ob; o++){
		if (SIGN * fita[o] < SIGN * fitb[o]){
			return false;
		}
	}
	for (int o = 0; o < ob; o++){
		if (SIGN * fita[o] > SIGN * fitb[o]){
			return true;
		}
	}
	return false;
}

void non_dominated_set(double **f, double **non_dom_set, int &ndsize, int M, int NPOP, int SIGN){
	int *same_solution = new int[NPOP];
	int *flag_dominated = new int[NPOP];
	int *same_group = new int[NPOP];
	for (int n = 0; n < NPOP; n++){
		same_solution[n] = 0;
		flag_dominated[n] = 0;
		same_group[n] = 0;
	}

	int same_group_count = 1;
	for (int b = 0; b < NPOP; b++){
		int flag = 0;
		for (int a = 0; a < NPOP; a++){
			if (a != b){
				if (dominate(f[a], f[b], M, SIGN)){
					flag_dominated[b] = 1;
				}
				if (same_solution[b] == 0 || same_solution[b] == same_group_count){
					int o = 0;
					for (; o < M; o++){
						if (f[a][o] != f[b][o]){
							break;
						}
					}
					if (o == M){
						same_solution[b] = same_group_count;
						same_solution[a] = same_group_count;
						flag = 1;
					}
				}
			}
		}
		if (flag == 1){
			same_group_count++;
		}
	}

	int nds = 0;
	for (int n = 0; n < NPOP; n++){
		if (flag_dominated[n] == 0 && same_solution[n] == 0){
			for (int o = 0; o < M; o++){
				non_dom_set[nds][o] = f[n][o];
			}
			nds++;

		}
		else if (flag_dominated[n] == 0 && same_solution[n] > 0){
			if (same_group[same_solution[n]] == 0){
				for (int o = 0; o < M; o++){
					non_dom_set[nds][o] = f[n][o];
				}
				same_group[same_solution[n]] = 1;
				nds++;
			}

		}
	}
	ndsize = nds;

	delete [] same_solution;
	delete [] flag_dominated;
	delete [] same_group;

}


int main(int argc, char *argv[]){
	std::chrono::system_clock::time_point  start, end;

	bool combination = false;
	
	if (strcmp(argv[1], "-num") != 0){
		cout << "file num " << endl;
		exit(1);
	}
	
	if (argv[2][0] == '-'){
		cout << "experimental number" << endl;
		exit(1);
	}
	init_genrand(atoi(argv[2]));


	int EPN = 6000;

	for (int i = 1; i < argc; i++){
		if (i % 2){
			if (argv[i][0] != '-'){
				cout << "format is disavailable" << endl;
				cout << "-... value -... value -... value ........" << endl;
				exit(1);
			}
		}
	}

	double knap_alpha = 0.0;
	char problem[50] = "wfg1";
	int ob = 10;
	int POPSIZE1 = 220;
	int MOEAD_H1 = 3;
	int T = POPSIZE1 / 10;

	double ALPHA_MAX = 1.0;
	double ALPHA_MIN = 1.0;
	double PENALTY = 5.0;
	double PENALTY2 = 50.0;
	int GEN = INT_MAX;
	int VALNUM = INT_MAX;

	double shita = M_PI / (2.0 * MOEAD_H1);
	double SBX_ETA = 20.0;
	double MUT_ETA = 20.0;
	double CO = 1.0;
	double mut = 1.0;
	int k_dtlz = 10;
	int k_fac = 2;
	int l_fac = 10;
	bool improve_lambda = false;
	bool normalize_flag = false;
	int POPSIZE2 = 0;
	int MOEAD_H2 = 1;
	char function[100] = "tchebycheff";

	bool ref_non_dom_flag = false;

	char max_min[10] = "min";

	int read_ind = 3;
	for (int read_ind = 0; read_ind < argc; read_ind++){
		if (strcmp(argv[read_ind], "-pro") == 0){
			if (strcmp(argv[read_ind + 1], "knapsack") == 0 ||
				strcmp(argv[read_ind + 1], "wfg1") == 0 || strcmp(argv[read_ind + 1], "wfg2") == 0 ||
				strcmp(argv[read_ind + 1], "wfg3") == 0 || strcmp(argv[read_ind + 1], "wfg4") == 0 ||
				strcmp(argv[read_ind + 1], "wfg5") == 0 || strcmp(argv[read_ind + 1], "wfg6") == 0 ||
				strcmp(argv[read_ind + 1], "wfg7") == 0 || strcmp(argv[read_ind + 1], "wfg8") == 0 ||
				strcmp(argv[read_ind + 1], "wfg9") == 0 ||
				strcmp(argv[read_ind + 1], "dtlz1") == 0 || strcmp(argv[read_ind + 1], "dtlz2") == 0 ||
				strcmp(argv[read_ind + 1], "dtlz3") == 0 || strcmp(argv[read_ind + 1], "dtlz4") == 0 ||
				strcmp(argv[read_ind + 1], "dtlz7") == 0 ||
				strcmp(argv[read_ind + 1], "maxwfg1") == 0 || strcmp(argv[read_ind + 1], "maxwfg2") == 0 ||
				strcmp(argv[read_ind + 1], "maxwfg3") == 0 || strcmp(argv[read_ind + 1], "maxwfg4") == 0 ||
				strcmp(argv[read_ind + 1], "maxwfg5") == 0 || strcmp(argv[read_ind + 1], "maxwfg6") == 0 ||
				strcmp(argv[read_ind + 1], "maxwfg7") == 0 || strcmp(argv[read_ind + 1], "maxwfg8") == 0 ||
				strcmp(argv[read_ind + 1], "maxwfg9") == 0 ||
				strcmp(argv[read_ind + 1], "maxdtlz1") == 0 || strcmp(argv[read_ind + 1], "maxdtlz2") == 0 ||
				strcmp(argv[read_ind + 1], "maxdtlz3") == 0 || strcmp(argv[read_ind + 1], "maxdtlz4") == 0 ||
				strcmp(argv[read_ind + 1], "maxdtlz7") == 0 ||
				strcmp(argv[read_ind + 1], "inverted_dtlz1") == 0 || strcmp(argv[read_ind + 1], "inverted_dtlz2") == 0 ||
				strcmp(argv[read_ind + 1], "zdt1") == 0 || strcmp(argv[read_ind + 1], "zdt2") == 0 ||
				strcmp(argv[read_ind + 1], "zdt3") == 0 || strcmp(argv[read_ind + 1], "zdt4") == 0 ||
				strcmp(argv[read_ind + 1], "mtest_dtlz") == 0 || strcmp(argv[read_ind + 1], "mtest_wfg") == 0 ||
				strcmp(argv[read_ind + 1], "test20160301") == 0 || strcmp(argv[read_ind + 1], "test20160303") == 0 ||
				strcmp(argv[read_ind + 1], "fonseca") == 0 || strcmp(argv[read_ind + 1], "dmp_twovar") == 0 ||
				strcmp(argv[read_ind + 1], "dmp") == 0 || strcmp(argv[read_ind + 1], "pdmp") == 0){
				strcpy(problem, argv[read_ind + 1]);
			}
			else{
				cout << "there is no problem." << endl;
				exit(1);
			}
		}
		else if (strcmp(argv[read_ind], "-obj") == 0){
			ob = atoi(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-pop") == 0){
			POPSIZE1 = atoi(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-dev") == 0){
			MOEAD_H1 = atoi(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-nei") == 0){
			T = atoi(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-fun") == 0){
			strcpy(function, argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-amax") == 0){
			ALPHA_MAX = atof(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-amin") == 0){
			ALPHA_MIN = atof(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-pena") == 0){
			PENALTY = atof(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-pena2") == 0){
			PENALTY2 = atof(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-shita") == 0){
			shita = atof(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-valnum") == 0){
			VALNUM = atoi(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-gen") == 0){
			GEN = atoi(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-sbxeta") == 0){
			SBX_ETA = atof(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-muteta") == 0){
			MUT_ETA = atof(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-co") == 0){
			CO = atof(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-mut") == 0){
			mut = atof(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-k_fac") == 0){
			k_fac = atoi(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-l_fac") == 0){
			l_fac = atoi(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-k_dtlz") == 0){
			k_dtlz = atoi(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-lam") == 0){
			if (strcmp(argv[read_ind + 1], "true") == 0){
				improve_lambda = true;
			}
		}
		else if (strcmp(argv[read_ind], "-maxmin") == 0){
			if (strcmp(argv[read_ind + 1], "max") == 0){
				strcpy(max_min, "max");
			}
		}
		else if (strcmp(argv[read_ind], "-refnondom") == 0){
			if (strcmp(argv[read_ind + 1], "true") == 0){
				ref_non_dom_flag = true;
			}
		}
		else if (strcmp(argv[read_ind], "-normalize") == 0){
			if (strcmp(argv[read_ind + 1], "true") == 0){
				normalize_flag = true;
			}
		}
		else if (strcmp(argv[read_ind], "-pop2") == 0){
			POPSIZE2 = atoi(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-dev2") == 0){
			MOEAD_H2 = atoi(argv[read_ind + 1]);
		}
		else if (strcmp(argv[read_ind], "-knapalpha") == 0){
			knap_alpha = atof(argv[read_ind + 1]);
		}

	}

	int point_num = 5;
	double **P = new double *[point_num];
	for (int p = 0; p < point_num; p++){
		P[p] = new double[ob];
	}

	
	double **extremepoint = new double *[ob];
	for (int i = 0; i < ob; i++){
		extremepoint[i] = new double[ob];
		for (int o = 0; o < ob; o++){
			extremepoint[i][o] = 0.0;
		}
	}
	
	double *interception = new double[ob];
	

	int k_val = k_fac;
	int l_val = l_fac;
	int item = k_val + l_val;


	if (strcmp("dtlz1", problem) == 0 || strcmp("dtlz2", problem) == 0 ||
		strcmp("dtlz3", problem) == 0 || strcmp("dtlz4", problem) == 0 || 
		strcmp("dtlz7", problem) == 0 ||
		strcmp("maxdtlz1", problem) == 0 || strcmp("maxdtlz2", problem) == 0 ||
		strcmp("maxdtlz3", problem) == 0 || strcmp("maxdtlz4", problem) == 0 ||
		strcmp("maxdtlz7", problem) == 0 ||
		strcmp("inverted_dtlz1", problem) == 0 || strcmp("inverted_dtlz2", problem) == 0){
		item = ob + k_dtlz - 1;
	}
	else if (strcmp("knapsack", problem) == 0){
		combination = true;

		char file_name[50];
		strcpy(file_name, problem);
		strcat(file_name, "_2_500");
		if (ob != 2){
			char ch[5];
			sprintf(ch, "%d", ob);
			strcat(file_name, "to");
			strcat(file_name, ch);
		}
		strcat(file_name, ".txt");
		knapsack_file_read(file_name, item, ob); //file_read
		sorting_profit_per_weight(item, ob); //q_j sorting
		repair_output(item, ob);
		int *cheeck = new int[item];
		for (int i = 0; i < item; i++){
			cheeck[i] = 1;
		}
		check_input_file_output(cheeck, item, ob);

	}
	else if (strcmp("zdt1", problem) == 0 || strcmp("zdt2", problem) == 0 ||
		strcmp("zdt3", problem) == 0 || strcmp("zdt4", problem) == 0){
		k_dtlz = 30;
		if (strcmp("zdt4", problem) == 0){
			k_dtlz = 10;
		}
		item = k_dtlz;
	}
	else if (strcmp("dmp_twovar", problem) == 0){
		k_dtlz = 2;
		item = 2;
	}
	else if (strcmp("dmp", problem) == 0 || strcmp("pdmp", problem) == 0){
		item = k_dtlz;
	}
	else if (strcmp("mtest_wfg", problem) == 0 || strcmp("mtest_dtlz", problem) == 0){
		point_num = 5;
		*P = new double[point_num];
		for (int i = 0; i < point_num; i++){
			P[i] = new double[ob];
		}
		ifstream fcp("cp8.txt");
		for (int p = 0; p < point_num; p++){
			for (int o = 0; o < ob; o++){
				fcp >> P[p][o];
			}
		}
		fcp.close();
		item = k_dtlz;
	}
	else if (strcmp("test20160301", problem) == 0){
		item = k_dtlz;
		combination = true;

	}
	else if (strcmp("test20160303", problem) == 0){
		item = k_dtlz;
		combination = true;

	}
	else if (strcmp("fonseca", problem) == 0){
		item = 2;
		combination = false;
	}
	else{
		*P = new double[point_num];
		for (int i = 0; i < point_num; i++){
			P[i] = new double[ob];
			for (int o = 0; o < ob; o++){
				P[i][o] = 0.0;
			}
		}

	}
	
	cout << "number: " << atoi(argv[2]) << endl;
	cout << "problem: " << problem << ", " << ob << "ob" << endl;
	cout << "popsize: " << POPSIZE1 << ", devision: " << MOEAD_H1 << ", neigborsize: " << T << endl;
	cout << "function: " << function << ", maxalpha: " << ALPHA_MAX << ", minalpha: " << ALPHA_MIN << endl;
	cout << "penalty1: " << PENALTY << ", penaly2: " << PENALTY2 << endl;
	cout << "shita: " << shita << endl;
	cout << "valnum: " << VALNUM << ", generation: " << GEN << endl << endl;
	cout << "sbx_eta: " << SBX_ETA << ", mut_eta: " << MUT_ETA << endl;
	cout << "crossover probabilty: " << CO << ", mutation probability: " << mut << " / " << item << endl;
	cout << "k_fac: " << k_fac << ", l_fac: " << l_fac << endl;
	cout << "k_dtlz: " << k_dtlz << endl;
	cout << "improve_lambda: " << improve_lambda << endl;
	
	double MUT = mut / item;

	char datafile[50];
	char outputfile[30] = "graph";
	char ccc[5] = "00";
	sprintf(ccc, "%d", atoi(argv[2]));
	strcat(outputfile, ccc);
	strcpy(datafile, outputfile);
	strcat(outputfile, ".txt");

	int N1 = POPSIZE1;
	int H1 = MOEAD_H1;

	int N2 = POPSIZE2;
	int H2 = MOEAD_H2;

	int N = N1 + N2;

	double **lambda_h1 = new double*[N1];
	for (int i = 0; i < N1; i++) lambda_h1[i] = new double[ob];

	

	/////////////////////////////////////////////////////////////


	//simplex-lattice design
	int *t_lambda = new int[ob];
	for (int i = 0; i < ob - 1; i++) t_lambda[i] = 0;

	
	bool end_flag = false;
	for (int i = 0; i < N1; i++){
		int temp = 0;
		for (int j = 0; j < ob - 1; j++) temp += t_lambda[j];
		t_lambda[ob - 1] = H1 - temp;

		bool assign_flag = true;
		for (int j = 0; j < ob; j++){
			assign_flag &= (t_lambda[ob - 1] >= 0);
		}

		if (assign_flag){
			for (int j = 0; j < ob; j++){
				lambda_h1[i][j] = (double)(t_lambda[j] / (double)H1);
			}

		}

		for (int j = ob - 1; j >= 1; j--){
			if (t_lambda[j] <= 0){ //continue affirmation
				if (j == 1) end_flag = true;
				int sum = 0;
				for (int k = 0; k < j; k++){
					sum += t_lambda[k];
					if (sum > H1){
						t_lambda[k] = 0;
					}
				}
			}
			else{
				t_lambda[j - 1] += 1;
				for (int k = j; k < ob - 1; k++){
					t_lambda[k] = 0;
				}
				break;
			}
		}
		if (end_flag) break;
	}




	double **lambda_h2 = new double*[N2];
	for (int i = 0; i < N2; i++) lambda_h2[i] = new double[ob];

	int *tt_lambda = new int[ob];
	for (int i = 0; i < ob - 1; i++) tt_lambda[i] = 0;


	bool end_flag2 = false;
	for (int i = 0; i < N2; i++){
		int temp = 0;
		for (int j = 0; j < ob - 1; j++) temp += tt_lambda[j];
		tt_lambda[ob - 1] = H2 - temp;

		bool assign_flag = true;
		for (int j = 0; j < ob; j++){
			assign_flag &= (tt_lambda[ob - 1] >= 0);
		}

		if (assign_flag){
			for (int j = 0; j < ob; j++){
				lambda_h2[i][j] = (double)(tt_lambda[j] / (double)H2);
			}

		}

		for (int j = ob - 1; j >= 1; j--){
			if (tt_lambda[j] <= 0){ //continue affirmation
				if (j == 1) end_flag2 = true;
				int sum = 0;
				for (int k = 0; k < j; k++){
					sum += tt_lambda[k];
					if (sum > H2){
						tt_lambda[k] = 0;
					}
				}
			}
			else{
				tt_lambda[j - 1] += 1;
				for (int k = j; k < ob - 1; k++){
					tt_lambda[k] = 0;
				}
				break;
			}
		}
		if (end_flag2) break;
	}

	for (int i = 0; i < N2; i++){
		for (int o = 0; o < ob; o++){
			lambda_h2[i][o] = ((lambda_h2[i][o] + (1.0 / (double)ob)) / 2.0);
		}
	}

	//if you use two layer this is deleted from here.
	double **lambda = new double *[N];
	for (int i = 0; i < N; i++) lambda[i] = new double[ob];

	int n12 = 0;
	for (n12 = 0; n12 < N1; n12++){
		for (int o = 0; o < ob; o++){
			lambda[n12][o] = lambda_h1[n12][o];
		}
	}
	for (; n12 < N; n12++){
		for (int o = 0; o < ob; o++){
			lambda[n12][o] = lambda_h2[n12 - N1][o];
		}
	}


	double **lambda2 = new double*[N];
	for (int i = 0; i < N; i++) lambda2[i] = new double[ob];




	for (int n = 0; n < N; n++){
		for (int o = 0; o < ob; o++){
			lambda2[n][o] = lambda[n][o];
		}
	}
	//to here.


	//if you use two layer this can be used as lambda.

	//double **lambda2 = new double*[N];
	//for (int i = 0; i < N; i++) lambda2[i] = new double[ob];


	/*for (int n = 0; n < N; n++){
		for (int o = 0; o < ob; o++){
			lambda[n][o] = 1.0 - lambda[n][o];
		}
	}

	for (int n = 0; n < N; n++){
		double norm = 0;
		for (int o = 0; o < ob; o++){
			norm += lambda[n][o];
		}
		for (int o = 0; o < ob; o++){
			lambda[n][o] = lambda[n][o] / norm;
		}
	}*/
	//for (int n = N / 2; n < N; n++){
	//	for (int o = 0; o < ob; o++){
	//		lambda[n][o] = lambda_h1[n - N / 2][o];
	//	}

	//}

	/*for (int n = 0; n < N; n++){
		for (int o = 0; o < ob; o++){
			lambda2[n][o] = lambda[n][o];
		}
	}*/
	//to here.


	

	for (int n = 0; n < N; n++){
		for (int o = 0; o < ob; o++){
			if (lambda2[n][o] == 0){
				lambda2[n][o] = 0.000001;
			}
		}
	}



	int **B = new int*[N];
	for (int i = 0; i < N; i++){
		B[i] = new int[T];
	}
	double *zmin = new double[ob];
	double *zmax = new double[ob];

	double *fmin = new double[ob];
	double *fmax = new double[ob];


	double **x = new double*[N];
	for (int n = 0; n < N; n++){
		x[n] = new double[item];
	}

	int **x_com = new int *[N];
	for (int n = 0; n < N; n++){
		x_com[n] = new int[item];
	}

	double **FV = new double*[N];
	for (int n = 0; n < N; n++){
		FV[n] = new double[ob];
	}
	double *y = new double[item];
	int *y_com = new int[item];
	double *y_fit = new double[ob];

	bool vector_chase_flag = false;
	if (atoi(argv[2]) == 0){
		vector_chase_flag = true;
	}

	double **non_dom_set = new double *[N];
	for (int i = 0; i < N; i++){
		non_dom_set[i] = new double[ob];
	}


	ofstream fout_ref_max;
	ofstream fout_ref_min;

	

	if (vector_chase_flag == true){
		fout_ref_max.open("ref_max0.txt");
		fout_ref_min.open("ref_min0.txt");

	}

	if (atoi(argv[2]) != 0){
		for (int i = 0; i < N; i++){
			int j = genrand_int31() % N;
			for (int o = 0; o < ob; o++){
				double t = lambda[i][o];
				lambda[i][o] = lambda[j][o];
				lambda[j][o] = t;
				t = lambda2[i][o];
				lambda2[i][o] = lambda2[j][o];
				lambda2[j][o] = t;
			}
		}

	}
	else{
		if (improve_lambda == false){
			ofstream vec("vec0.txt");
			for (int i = 0; i < N; i++){
				for (int o = 0; o < ob; o++){
					vec << lambda[i][o] << "\t";
				}
				vec << endl;
			}
		}
		else{
			ofstream vec("vec0.txt");
			for (int i = 0; i < N; i++){
				for (int o = 0; o < ob; o++){
					vec << lambda2[i][o] << "\t";
				}
				vec << endl;
			}
		}
	}

	cout << N << " " << T << endl;
	for (int o = 0; o < ob; o++){
		cout << lambda[N - 1][o] << " ";
	}
	cout << endl;

	///************************************Step1 Initialization***********************************/
	//Step1.1

	//Step1.2

	t_neighborhood(B, lambda, N, T, ob);

	//Step1.3

	if (combination == false){
		for (int n = 0; n < N; n++){
			for (int i = 0; i < item; i++){
				x[n][i] = nextDoubleIE();
			}
		}
	}
	else{

		if (strcmp("test20160301", problem) == 0 || strcmp("test20160303", problem) == 0){
			for (int n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					x_com[n][i] = genrand_int31() % ob + 1;
				}
			}

		}
		else if (strcmp("test20160303", problem) == 0){
			for (int n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					x_com[n][i] = genrand_int31() % k_fac + 1;
				}
			}

		}
		else{
			for (int n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					x_com[n][i] = genrand_int31() % 2;
				}
			}
		}
	}


	if (combination == false){
		for (int n = 0; n < N; n++){
			fitness(x[n], FV[n], problem, ob, k_val, l_val, k_dtlz, P, point_num);
		
		}
	}
	else{
		if (strcmp(problem, "knapsack") == 0){
			if (strcmp(max_min, "min") == 0){
				for (int n = 0; n < N; n++){
					minus_fitness_knap(x_com[n], FV[n], item, ob, knap_alpha);
				}
			}
			else{
				for (int n = 0; n < N; n++){
					fitness_knap(x_com[n], FV[n], item, ob, knap_alpha);
				}
			}
		}
		else if (strcmp("test20160301", problem) == 0){
			for (int n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					test20160301(x_com[n], FV[n], ob, item);
				}
			}
		}
		else if (strcmp("test20160303", problem) == 0){
			for (int n = 0; n < N; n++){
				for (int i = 0; i < item; i++){
					test20160303(x_com[n], FV[n], ob, item, k_fac);
				}
			}
		}
		else{
			exit(1);
		}
	}

	
	ofstream fffout("dtlz-ini.dat");

	for (int i = 0; i < N; i++){
	for (int o = 0; o < ob; o++){
	fffout << FV[i][o] << "\t";
	}
	fffout << endl;

	}

	fffout.close();

	
	

	//Step1.4
	for (int o = 0; o < ob; o++){
		zmin[o] = 100000.0;
		fmin[o] = 1000000.0;
	}
	for (int o = 0; o < ob; o++){
		zmax[o] = 0.0;
		fmax[o] = 0.0;
	}
	for (int n = 0; n < N; n++){
		max_z(FV[n], zmax, ob, ALPHA_MAX, fmax);
	}

	for (int n = 0; n < N; n++){
		min_z(FV[n], zmin, ob, ALPHA_MIN, fmin);
	}
	int plus_minus = -1;

	if (ref_non_dom_flag == true){
		int ndsize = 0;
		int SIGN = 1;
		if (strcmp(max_min, "min") == 0){
			SIGN = -1;
		}
		non_dominated_set(FV, non_dom_set, ndsize, ob, N, SIGN);
		for (int n = 0; n < ndsize; n++){
			max_z(non_dom_set[n], zmax, ob, ALPHA_MAX, fmax);
		}

		for (int n = 0; n < ndsize; n++){
			min_z(non_dom_set[n], zmin, ob, ALPHA_MIN, fmin);
		}
	}


	/************************************Step2 Update***********************************/

	//progress_output
	char detaildata[100] = "./";
	strcat(detaildata, datafile);
	strcat(detaildata, "/");
	strcat(detaildata, "eva");
	char gennumch[8] = "00";
	sprintf(gennumch, "%d", 0);
	strcat(detaildata, gennumch);
	strcat(detaildata, ".txt");


	double *normalized_FV = new double[ob];
	double *normalized_y_fit = new double[ob];
	double *normalized_zmin = new double[ob];
	double *normalized_zmax = new double[ob];




	double alpha = ALPHA_MAX;
	double alpha_max = ALPHA_MAX;
	double alpha_min = ALPHA_MIN;
	int divide = VALNUM / N;
	//cout << divide << endl;
	int finiflag = 0;


	//output
	//from here
	//char range_maxminoutputfile[100][40];
	//char tempchar2[5] = "00";
	//for (int o = 0; o < ob; o++){
	//	strcpy(range_maxminoutputfile[o], "range_max_min");
	//	sprintf(tempchar2, "%d", atoi(argv[2]));
	//	strcat(range_maxminoutputfile[o], "_graph");
	//	strcat(range_maxminoutputfile[o], tempchar2);
	//	strcat(range_maxminoutputfile[o], "_");
	//	sprintf(tempchar2, "%d", o + 1);
	//	strcat(range_maxminoutputfile[o], tempchar2);
	//	strcat(range_maxminoutputfile[o], "ob.txt");
	//}


	//char range_nondom_maxminoutputfile[100][40];
	//for (int o = 0; o < ob; o++){
	//	strcpy(range_nondom_maxminoutputfile[o], "range_nondom_max_min");
	//	sprintf(tempchar2, "%d", atoi(argv[2]));
	//	strcat(range_nondom_maxminoutputfile[o], "_graph");
	//	strcat(range_nondom_maxminoutputfile[o], tempchar2);
	//	strcat(range_nondom_maxminoutputfile[o], "_");
	//	sprintf(tempchar2, "%d", o + 1);
	//	strcat(range_nondom_maxminoutputfile[o], tempchar2);
	//	strcat(range_nondom_maxminoutputfile[o], "ob.txt");
	//}


	//output
	//to here

	for (int g = 0; g < GEN; g++){

		//only minimization
		int ndsize = 0;
		int SIGN = -1;

		non_dominated_set(FV, non_dom_set, ndsize, ob, N, SIGN);

		//output
		//from here
	/*	char genoutputfile[30] = "gen";
		char genoutputfile_objective[100][40];
		char tempchar[5] = "00";
		sprintf(tempchar, "%d", g);
		strcat(genoutputfile, tempchar);
		strcat(genoutputfile, "_graph");
		sprintf(tempchar, "%d", atoi(argv[2]));
		strcat(genoutputfile, tempchar);
		for (int o = 0; o < ob; o++){
			strcpy(genoutputfile_objective[o], genoutputfile);
			sprintf(tempchar, "%d", o + 1);
			strcat(genoutputfile_objective[o], "_");
			strcat(genoutputfile_objective[o], tempchar);
			strcat(genoutputfile_objective[o], "ob.txt");
		}

		strcat(genoutputfile, ".txt");
		ofstream genoutput(genoutputfile);
		
		genoutput << setprecision(20);
		for (int nd = 0; nd < ndsize; nd++){
			for (int o = 0; o < ob - 1; o++){
				genoutput << non_dom_set[nd][o] << "\t";
			}
			genoutput << non_dom_set[nd][ob - 1];
			genoutput << endl;
		}
		genoutput.close();
		double *nondom_min = new double[ob];
		double *nondom_max = new double[ob];
		for (int o = 0; o < ob; o++){
			nondom_min[o] = DBL_MAX;
			nondom_max[o] = DBL_MIN;
			ofstream genoutputob(genoutputfile_objective[o]);
			ofstream range_nondom_maxminoutput(range_nondom_maxminoutputfile[o], ios::app);
			ofstream range_maxminoutput(range_maxminoutputfile[o], ios::app);
			range_nondom_maxminoutput << setprecision(20);
			range_maxminoutput << setprecision(20);
			genoutputob << setprecision(20);
			range_maxminoutput << setprecision(20);
			for (int nd = 0; nd < ndsize; nd++){
				genoutputob << non_dom_set[nd][o] << endl;
				if (non_dom_set[nd][o] > nondom_max[o]){
					nondom_max[o] = non_dom_set[nd][o];
				}
				if (non_dom_set[nd][o] < nondom_min[o]){
					nondom_min[o] = non_dom_set[nd][o];
				}
			}
			range_nondom_maxminoutput <<g << "\t" <<  (nondom_max[o] - nondom_min[o]) << endl;
			range_nondom_maxminoutput.close();
			range_maxminoutput << g << "\t" << (nondom_max[o] - zmin[o]) << endl;
			range_maxminoutput.close();
			genoutputob.close();
		}
		delete[] nondom_min;
		delete[] nondom_max;
*/
		//output
		//to here

		double *realmin = new double[ob];
		double *realmax = new double[ob];
		for (int o = 0; o < ob; o++){
			realmin[o] = zmin[o];
			realmax[o] = non_dom_set[0][o];
		}
		for (int o = 0; o < ob; o++){
			zmax[o] = 0.0;
			fmax[o] = 0.0;
		}
		for (int n = 0; n < N; n++){
			max_z(FV[n], zmax, ob, ALPHA_MAX, fmax);
		}
		for (int ns = 0; ns < ndsize; ns++){
			for (int o = 0; o < ob; o++){
				//if (realmax[o] < non_dom_set[ns][o]){
					realmax[o] = zmax[o];
				//}

			}
		}
		//only minimization


		

		CalcExtremePoint(extremepoint, non_dom_set, realmax, realmin, ob, N);

		CalcInterception(ob, extremepoint, realmin, realmax, interception);



	
	/*	CalcExtremePoint(extremepoint, FV, zmax, zmin, ob, N);

		CalcInterception(ob, extremepoint, zmin, zmax, interception);*/

		delete[] realmin;
		delete[] realmax;
		for (int i = 0; i < N; i++){
			if (g * N + i < VALNUM){
			
				
			
				//Step2.1 Reproduction
				int k = genrand_int31() % T;
				int l = genrand_int31() % T;
				
				
				if (combination == false){
					make_new_pop(x[B[i][k]], x[B[i][l]], y, item, MUT, CO, MUT_ETA, SBX_ETA);
				}
				else{
					if (strcmp(problem, "knapsack") == 0){
						make_new_pop(x_com[B[i][k]], x_com[B[i][l]], y_com, item, MUT, CO);
					}
					else if (strcmp(problem, "test20160301") == 0){
						make_new_pop_20160301(x_com[B[i][k]], x_com[B[i][l]], y_com, item, MUT, CO, ob);
					}
					else if (strcmp(problem, "test20160303") == 0){
						make_new_pop_20160301(x_com[B[i][k]], x_com[B[i][l]], y_com, item, MUT, CO, k_fac);
					}
					else{
						exit(1);
					}
				}
			
				//Step2.2 Improvement
				if (combination == false){
					fitness(y, y_fit, problem, ob, k_val, l_val, k_dtlz, P, point_num);
			
				}
				else{
					if (strcmp(problem, "knapsack") == 0){
						if (strcmp(max_min, "min") == 0){
							minus_fitness_knap(y_com, y_fit, item, ob, knap_alpha);
							
						}
						else{
							fitness_knap(y_com, y_fit, item, ob, knap_alpha);
						}
					}
					else if (strcmp(problem, "test20160301") == 0){
						test20160301(y_com, y_fit, ob, item);
					}
					else if (strcmp(problem, "test20160303") == 0){
						test20160303(y_com, y_fit, ob, item, k_fac);
					}
					else{
						exit(1);
					}
					
				}


				//Step2.3 Update of z
				if (strcmp(max_min, "min") == 0){
					min_z(y_fit, zmin, ob, ALPHA_MIN, fmin);
					
					for (int o = 0; o < ob; o++){
						zmax[o] = 0.0;
						fmax[o] = 0.0;
					}
					for (int n = 0; n < N; n++){
						max_z(FV[n], zmax, ob, ALPHA_MAX, fmax);
					}

				//	max_z(y_fit, zmax, ob, ALPHA_MAX, fmax);

				}
				else{
					max_z(y_fit, zmax, ob, ALPHA_MAX, fmax);
					
					for (int o = 0; o < ob; o++){
						zmin[o] = 100000.0;
						fmin[o] = 100000.0;
					}

					for (int n = 0; n < N; n++){
						min_z(FV[n], zmin, ob, ALPHA_MIN, fmin);
					}
				//	min_z(y_fit, zmin, ob, ALPHA_MIN, fmin);

				}

				//Step2.4 Update of Neighboring Solutions
			

			
			
			
				for (int j = 0; j < T; j++){
					int jejeje = B[i][j];
				
					for (int o = 0; o < ob; o++){
						normalized_FV[o] = (FV[jejeje][o] - zmin[o]) / (interception[o] - zmin[o]);
						normalized_y_fit[o] = (y_fit[o] - zmin[o]) / (interception[o] - zmin[o]);
						normalized_zmin[o] = 0.0;
						normalized_zmax[o] = 1.0;					
					}
					if (strcmp(function, "normalized_pbi") == 0){
						for (int o = 0; o < ob; o++){
							lambda2[jejeje][o] = interception[o] - zmin[o];
						}
					}

				
					double shita_temp = shita;
					if (strcmp(function, "fix_double_pbi") == 0 || strcmp(function, "fix_double_inverted_pbi") == 0
						|| strcmp(function, "quadratic_pbi") == 0 || strcmp(function, "quadratic_inverted_pbi") == 0){
						double diagonal = 0.0;
						for (int o = 0; o < ob; o++){
							diagonal += (fmax[o] - fmin[o]);
						}
						shita = /*(1.0 + 100.0 * (1.0 - ((double)g + 1.0) / ((double)(VALNUM / N) + 1.0))) */ shita_temp * diagonal / (double)ob / (double)H1;
						/*if (!(i % N) && !(j % T)){
							cout << (VALNUM / N) << " " <<(1.0 + 100.0 * (1.0 - ((double)g + 1.0) / ((double)(VALNUM / N) + 1.0))) << endl;
						}*/
					}
					if (function_valuation(problem, function,
						FV[jejeje], y_fit, 
						improve_lambda,lambda[jejeje],lambda2[jejeje], 
						zmax, zmin, 
						ob, 
						PENALTY, PENALTY2, shita, max_min) && (!normalize_flag)){
					
						if (combination == false){
							for (int it = 0; it < item; it++){
								x[jejeje][it] = y[it];
							}
						}
						else{
							for (int it = 0; it < item; it++){
								x_com[jejeje][it] = y_com[it];
							}
						}
						for (int o = 0; o < ob; o++){
							FV[jejeje][o] = y_fit[o];
						}
						
					}
					else if (function_valuation(problem, function,
						normalized_FV, normalized_y_fit,
						improve_lambda, lambda[jejeje], lambda2[jejeje],
						normalized_zmax, normalized_zmin,
						ob,
						PENALTY, PENALTY2, shita, max_min) && normalize_flag){
						
						if (combination == false){
							for (int it = 0; it < item; it++){
								x[jejeje][it] = y[it];
							}
						}
						else{
							for (int it = 0; it < item; it++){
								x_com[jejeje][it] = y_com[it];
							}
						}
						for (int o = 0; o < ob; o++){
							FV[jejeje][o] = y_fit[o];
						}
					}
			

					shita = shita_temp;
				}
			
				if (vector_chase_flag == true){
					if (!((g * N + i) % 105)){
						for (int o = 0; o < ob - 1; o++){
							fout_ref_max << zmax[o] << "\t";
						}
						fout_ref_max << zmax[ob - 1];
						fout_ref_max << endl;
						for (int o = 0; o < ob - 1; o++){
							fout_ref_min << zmin[o] << "\t";
						}
						fout_ref_min << zmin[ob - 1];
						fout_ref_min << endl;
					}
				}
				if (!((g * N + i) % 100000)){
					cout << g * N + i << endl;
					
				}

				//progress_output
				if (!((g * N + i) % 100)){
					strcpy(detaildata, "./");
					strcat(detaildata, datafile);
					strcat(detaildata, "/");
					strcat(detaildata, "eva");

					sprintf(gennumch, "%d", g * N + i);
					strcat(detaildata, gennumch);
					strcat(detaildata, ".txt");


				}
			}
			else{

				finiflag = 1;
				break;
			}

		}

		if (finiflag == 1){
			break;
		}
	}





	ofstream fout(outputfile);
	
	if (vector_chase_flag = true){
		for (int o = 0; o < ob - 1; o++){
			fout_ref_max << zmax[o] << "\t";
		}
		fout_ref_max << zmax[ob - 1];
		fout_ref_max << endl;
		for (int o = 0; o < ob - 1; o++){
			fout_ref_min << zmin[o] << "\t";
		}
		fout_ref_min << zmin[ob - 1];
		fout_ref_min << endl;
		fout_ref_max.close();
		fout_ref_min.close();
	}
	fout << setprecision(20);
	for (int i = 0; i < N; i++){
		for (int o = 0; o < ob - 1; o++){
			fout << FV[i][o] << "\t";
		}
		fout << FV[i][ob - 1];
		fout << endl;
	}

	fout.close();



	if (strcmp("dmp_twovar", problem) == 0){
		char varoutfile[30] = "var";
		char cccc[5] = "00";
		sprintf(cccc, "%d", atoi(argv[2]));
		strcat(varoutfile, cccc);
		strcat(varoutfile, ".txt");
		ofstream varout(varoutfile);
		for (int i = 0; i < N; i++){
			for (int it = 0; it < item - 1; it++){
				varout << 100.0 * x[i][it] - 50.0 << "\t";
			}
			varout << 100.0 * x[i][item - 1] - 50.0;
			varout << endl;
		}
		varout.close();
	}


	if (strcmp("dmp", problem) == 0 || strcmp("pdmp", problem) == 0){
		char varoutfile[30] = "var";
		char cccc[5] = "00";
		sprintf(cccc, "%d", atoi(argv[2]));
		strcat(varoutfile, cccc);
		strcat(varoutfile, ".txt");
		ofstream varout(varoutfile);
		for (int i = 0; i < N; i++){
			double x_v1 = 0.0;
			double v1norm = 0.0;
			double x_v2 = 0.0;
			double v2norm = 0.0;
			for (int it = 0; it < item; it++){
				if (!(it % 2)){
					x_v1 += 100.0 * x[i][it] - 50.0;
					v1norm += 1.0;
				}
				else{
					x_v2 += 100.0 * x[i][it] - 50.0;
					v2norm += 1.0;
				}
			}
			varout << x_v1 / v1norm << "\t" << x_v2 / v2norm;
			varout << endl;
		}
		varout.close();
	}



	delete[] normalized_FV;
	delete[] normalized_y_fit;
	delete[] normalized_zmin;
	delete[] normalized_zmax;
	
	//foutdetail3.close();

	for (int i = 0; i < ob; i++){
		delete[] extremepoint[i];
	}
	delete[] extremepoint;

	delete[] interception;

	return 0;
}

