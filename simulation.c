#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

const gsl_rng *gBaseRand;

int main(int argc, char* argv[]) {

	double s1, s2, p_asy, p_s, pref, rpref_S, rpref_A;            //Model parameters
	double avg_h_cumu = 0., curr_h_cumu = 0.;                     //Recording HL
	double avg_gen = 0., curr_gen = 0., avg_fix_gen = 0.;         //Recording generations
	double new_VS = 0., new_VA = 0., rfit_VS = 0., rfit_VA = 0.;  //Used for updating
	double subpop_freq_VS = 0., subpop_freq_VA = 0.;              //Frequencies of V in subpop

	int SEED, N, TIME_STEPS, REPS;       //Simulation setups
	int total_S, total_A;
	int poly_case = 0, curr_poly = 0;    //Recording polymorphism
	int fix_case = 0, curr_fix = 0;      //Recording fixation
	int curr_rep = 0, t = 0;
	int first_host = 0, sampled_VS = 0, sampled_VA = 0; //Used for updating
	int continue_run = 1;


	N = atof(argv[1]);
	s1 = atof(argv[2]);
	s2 = atof(argv[3]);
	p_asy = atof(argv[4]);
	pref = atof(argv[5]);
	REPS = atof(argv[6]);
	SEED = atof(argv[7]);
	printf("running rewrite\n");
	printf("N = %d\n", N);
	printf("s1 = %f\n", s1);
	printf("s2 = %f\n", s2);
	printf("p_asy = %f\n", p_asy);
	printf("pref = %f\n", pref);
	printf("reps = %d\n", REPS);

	TIME_STEPS = 10000000;
	//TIME_STEPS = 1000;      //For testing

	p_s = 1. - p_asy;
	rpref_S = pref; //Absolute preference
	rpref_A = pref;

	total_S = round(p_s * N);
	total_A = round(p_asy * N);

	clock_t begin = clock();
	gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(gBaseRand, SEED); 

	
	//LOOP: SIMULATE REPS
	while (curr_rep < REPS) {
		//printf("curr rep: %d\n", curr_rep);
		//Reset results and flags for one-run
		curr_rep += 1;

		curr_h_cumu = curr_gen = 0.;
		curr_fix = t = 0;
		continue_run = 1;

		//Reset frequencies; host of first V introduced is sampled randomly
		first_host = gsl_ran_binomial(gBaseRand, p_asy, 1);
		if (first_host == 0) {
			subpop_freq_VS = 1. / total_S; //first on S
			subpop_freq_VA = 0.;
		} else {
			subpop_freq_VA = 1. / total_A; //first on A
			subpop_freq_VS = 0.;
		}


		//LOOP: SIMULATE ONE REP
		while ((t < TIME_STEPS) && (continue_run == 1)) {
			t += 1;


			//Preference
			new_VS = subpop_freq_VS * rpref_S + subpop_freq_VA * (1. - rpref_S);
			new_VA = subpop_freq_VA * rpref_A + subpop_freq_VS * (1. - rpref_A);
			subpop_freq_VS = new_VS;
			subpop_freq_VA = new_VA;


			//Selection
			rfit_VS = 1. / (1. + (s2 * (1. - subpop_freq_VS)));
			rfit_VA = (1. + s1) / (1. + (s1 * subpop_freq_VA));
			new_VS = rfit_VS * subpop_freq_VS;
			new_VA = rfit_VA * subpop_freq_VA;
			subpop_freq_VS = new_VS;
			subpop_freq_VA = new_VA;


			//Sampling
			sampled_VS = gsl_ran_binomial(gBaseRand, subpop_freq_VS, total_S);
			sampled_VA = gsl_ran_binomial(gBaseRand, subpop_freq_VA, total_A);
			if (total_S == 0) {
				subpop_freq_VS = 0.;
			} else {
				subpop_freq_VS = (double) sampled_VS / total_S;
			}
			if (total_A == 0) {
				subpop_freq_VA = 0.;
			} else {
				subpop_freq_VA = (double) sampled_VA / total_A;
			}


			//Update subdivided heterozygosity
			curr_h_cumu += (subpop_freq_VS * (1. - subpop_freq_VS) * p_s);
			curr_h_cumu += (subpop_freq_VS * (1. - subpop_freq_VA) * p_asy);


			//Set flags
			if ((subpop_freq_VS == 0.) && (subpop_freq_VA == 0.)) {
				continue_run = 0; //V extinct
			} else if ((subpop_freq_VS == 1.) && (subpop_freq_VA == 1.)) {
				continue_run = 0; //V fix
				curr_fix = 1;
			}
		}
		//END OF LOOP: SIMULATE ONE REP

		//Record results of one run
		poly_case += continue_run; //Stable polymorphism if V not extinct/fix after T
		avg_h_cumu += curr_h_cumu;
		//printf("%f\n", avg_h_cumu);
		//printf("%f\n", curr_h_cumu);
		avg_gen += curr_gen;
		fix_case += curr_fix;
		if (curr_fix == 1) {
			avg_fix_gen += curr_gen;
		}

	}
	//END OF LOOP: SIMULATE REPS


	printf("avg_gen = %f\n", avg_gen / REPS);
	printf("avg_h = %f\n", avg_h_cumu / REPS);
	printf("prob(poly) = %f\n", (double)poly_case / REPS);
	printf("prob(fix) = %f\n", (double)fix_case / REPS);
	printf("avg_fix_gen = %f\n", avg_fix_gen / fix_case);

	clock_t end = clock();
	double total_time = (double) (end - begin) / CLOCKS_PER_SEC;
	printf("seconds took: %f\n", total_time);

	return 0;


}