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

	double s1, s2, p_asy1, p_asy2, p_s1, p_s2, pref1, pref2;      //Model parameters
	double avg_h_cumu = 0., curr_h_cumu = 0.;                     //Recording HL
	double avg_gen = 0., curr_gen = 0., avg_fix_gen = 0.;         //Recording generations
	double new_VS = 0., new_VA = 0., rfit_VS = 0., rfit_VA = 0.;  //Used for updating
	double subpop_freq_VS1 = 0., subpop_freq_VA1 = 0.;            //Frequencies of V in subpop
	double subpop_freq_VS2 = 0., subpop_freq_VA2 = 0.;            //Frequencies of V in subpop
	double fV1 = 0., fV2 = 0.;                                    //Frequencies of V in metapop
	double rpref_S1, rpref_S2, rpref_A1, rpref_A2;                //For relative preferences
	double migration = 0., no_migration = 0.;

	int SEED, N, TIME_STEPS, REPS;       //Simulation setups
	int total_S1, total_A1, total_S2, total_A2;
	int poly_case = 0, curr_poly = 0;    //Recording polymorphism
	int fix_case = 0, curr_fix = 0;      //Recording fixation
	int curr_rep = 0, t = 0;
	int first_host = 0, first_deme = 0;  //For sampling first V
	int sampled_VS = 0, sampled_VA = 0;  //Used for updating
	int continue_run = 1;


	N = atof(argv[1]);
	s1 = atof(argv[2]);
	s2 = atof(argv[3]);
	p_asy1 = atof(argv[4]); // needs to be 0.5
	p_asy2 = atof(argv[5]); // needs to be 0.5
	pref1 = atof(argv[6]);
	pref2 = atof(argv[7]);
	migration = atof(argv[8]);
	REPS = atof(argv[9]);
	SEED = atof(argv[10]);
	printf("running rewrite\n");
	printf("N = %d\n", N);
	printf("s1 = %f\n", s1);
	printf("s2 = %f\n", s2);
	printf("p_asyA = %f\n", p_asy1);
	printf("p_asyB = %f\n", p_asy2);
	printf("prefA = %f\n", pref1);
	printf("prefB = %f\n", pref2);
	printf("mig = %f\n", migration);
	printf("reps = %d\n", REPS);

	TIME_STEPS = 10000000;
	//TIME_STEPS = 1000;      //For testing

	p_s1 = 1. - p_asy1;
	p_s2 = 1. - p_asy2;
	rpref_S1 = 2.0 * pref1 * p_s1; //relative preference: to correct for splitted to two demes and random interactions
	rpref_S2 = 2.0 * pref2 * p_s2;
	rpref_A1 = 2.0 * pref1 * p_asy1;
	rpref_A2 = 2.0 * pref2 * p_asy2;

	total_S1 = round(p_s1 * N);
	total_A1 = round(p_asy1 * N);
	total_S2 = round(p_s2 * N);
	total_A2 = round(p_asy2 * N);

	no_migration = 1. - migration;

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

		//Reset frequencies; host and deme of first V introduced is sampled randomly
		first_deme = gsl_ran_binomial(gBaseRand, 0.5, 1);

		if (first_deme == 0) {
			//Introduce in deme1
			subpop_freq_VS2 = 0.;
			subpop_freq_VA2 = 0.;
			first_host = gsl_ran_binomial(gBaseRand, p_asy1, 1);
			if (first_host == 0) {
				subpop_freq_VS1 = 1. / total_S1; //first on S
				subpop_freq_VA1 = 0.;
			} else {
				subpop_freq_VS1 = 0.;
				subpop_freq_VA1 = 1. / total_A1; //first on A
			}
		} else {
			//introduce in deme2
			subpop_freq_VS1 = 0.;
			subpop_freq_VA1 = 0.;
			first_host = gsl_ran_binomial(gBaseRand, p_asy2, 1);
			if (first_host == 0) {
				subpop_freq_VS2 = 1. / total_S2; //first on S
				subpop_freq_VA2 = 0.;
			} else {
				subpop_freq_VS2 = 0.;
				subpop_freq_VA2 = 1. / total_A2; //first on A
			}
		}

		//printf("init, %.8f, %.8f, %.8f, %.8f\n", subpop_freq_VS1, subpop_freq_VA1, subpop_freq_VS2, subpop_freq_VA2);


		//LOOP: SIMULATE ONE REP
		while ((t < TIME_STEPS) && (continue_run == 1)) {
			t += 1;

			//Preference, within each deme
			new_VS = subpop_freq_VS1 * rpref_S1 + subpop_freq_VA1 * (1. - rpref_S1);
			new_VA = subpop_freq_VA1 * rpref_A1 + subpop_freq_VS1 * (1. - rpref_A1);
			subpop_freq_VS1 = new_VS;
			subpop_freq_VA1 = new_VA;

			new_VS = subpop_freq_VS2 * rpref_S2 + subpop_freq_VA2 * (1. - rpref_S2);
			new_VA = subpop_freq_VA2 * rpref_A2 + subpop_freq_VS2 * (1. - rpref_A2);
			subpop_freq_VS2 = new_VS;
			subpop_freq_VA2 = new_VA;

			//printf("pref, %.8f, %.8f, %.8f, %.8f\n", subpop_freq_VS1, subpop_freq_VA1, subpop_freq_VS2, subpop_freq_VA2);


			//Selection, within each deme
			rfit_VS = 1. / (1. + (s2 * (1. - subpop_freq_VS1)));
			rfit_VA = (1. + s1) / (1. + (s1 * subpop_freq_VA1));
			new_VS = rfit_VS * subpop_freq_VS1;
			new_VA = rfit_VA * subpop_freq_VA1;
			subpop_freq_VS1 = new_VS;
			subpop_freq_VA1 = new_VA;
			
			rfit_VS = 1. / (1. + (s2 * (1. - subpop_freq_VS2)));
			rfit_VA = (1. + s1) / (1. + (s1 * subpop_freq_VA2));
			new_VS = rfit_VS * subpop_freq_VS2;
			new_VA = rfit_VA * subpop_freq_VA2;
			subpop_freq_VS2 = new_VS;
			subpop_freq_VA2 = new_VA;

			//printf("sele, %.8f, %.8f, %.8f, %.8f\n", subpop_freq_VS1, subpop_freq_VA1, subpop_freq_VS2, subpop_freq_VA2);


			//Migration between demes
			fV1 = (subpop_freq_VA1 * p_asy1) + (subpop_freq_VS1 * p_s1);
			fV2 = (subpop_freq_VA2 * p_asy2) + (subpop_freq_VS2 * p_s2);
			subpop_freq_VS1 = ((subpop_freq_VS1 * p_s1 * no_migration) + (fV2 * migration * p_s1)) / p_s1;
			subpop_freq_VA1 = ((subpop_freq_VA1 * p_asy1 * no_migration) + (fV2 * migration * p_asy1)) / p_asy1;
			subpop_freq_VS2 = ((subpop_freq_VS2 * p_s2 * no_migration) + (fV1 * migration * p_s2)) / p_s2;
			subpop_freq_VA2 = ((subpop_freq_VA2 * p_asy2 * no_migration) + (fV1 * migration * p_asy2)) / p_asy2;

			//printf("migr, %.8f, %.8f, %.8f, %.8f\n", subpop_freq_VS1, subpop_freq_VA1, subpop_freq_VS2, subpop_freq_VA2);


			//Sampling
			sampled_VS = gsl_ran_binomial(gBaseRand, subpop_freq_VS1, total_S1);
			sampled_VA = gsl_ran_binomial(gBaseRand, subpop_freq_VA1, total_A1);
			subpop_freq_VS1 = (double) sampled_VS / total_S1;
			subpop_freq_VA1 = (double) sampled_VA / total_A1;

			sampled_VS = gsl_ran_binomial(gBaseRand, subpop_freq_VS2, total_S2);
			sampled_VA = gsl_ran_binomial(gBaseRand, subpop_freq_VA2, total_A2);
			subpop_freq_VS2 = (double) sampled_VS / total_S2;
			subpop_freq_VA2 = (double) sampled_VA / total_A2;	

			//printf("samp, %.8f, %.8f, %.8f, %.8f\n", subpop_freq_VS1, subpop_freq_VA1, subpop_freq_VS2, subpop_freq_VA2);		


			//Update subdivided heterozygosity
			curr_h_cumu += (subpop_freq_VS1 * (1. - subpop_freq_VS1) * p_s1 * 0.5);
			curr_h_cumu += (subpop_freq_VA1 * (1. - subpop_freq_VA1) * p_asy1 * 0.5);
			curr_h_cumu += (subpop_freq_VS2 * (1. - subpop_freq_VS2) * p_s2 * 0.5);
			curr_h_cumu += (subpop_freq_VA2 * (1. - subpop_freq_VA2) * p_asy2 * 0.5);


			//Set flags
			if ((subpop_freq_VS1 == 0.) && (subpop_freq_VA1 == 0.) && 
					(subpop_freq_VS2 == 0.) && (subpop_freq_VA2 == 0.)) {
				continue_run = 0; //V extinct
			} else if ((subpop_freq_VS1 == 1.) && (subpop_freq_VA1 == 1.) &&
					(subpop_freq_VS2 == 1.) && (subpop_freq_VA2 == 1.)) {
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