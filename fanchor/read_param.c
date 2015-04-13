#include "lb.h"

//Read in parameters from input FILE
void read_param()
{
	
	FILE *param;
	param = fopen("param.in", "r");
	fscanf(param, "newrun %d\n", &newrun_on);
	fscanf(param, "simul_evol %d\n", &simul_evol_on);
	
	fscanf(param, "particle_on %d\n", &particle_on);	
	fscanf(param, "wall_on %d\n", &wall_on);
	fscanf(param, "wall_anchoring_on %d\n", &wall_anchoring_on);
	fscanf(param, "W_bot %lf\n", &W_bot);
	fscanf(param, "W_top %lf\n", &W_top);
	fscanf(param, "nbot %lf %lf %lf\n", &nbot[0], &nbot[1], &nbot[2]);
	fscanf(param, "ntop %lf %lf %lf\n", &ntop[0], &ntop[1], &ntop[2]);
	
	fscanf(param, "flow_on %d\n", &flow_on);
	fscanf(param, "Q_on %d\n", &Q_on);
	fscanf(param, "n_evol_Q %d\n", &n_evol_Q);
	
	fscanf(param, "t_max %d\n", &t_max);
	fscanf(param, "t_print %d\n", &t_print);
	
	fscanf(param, "tau_f %le\n", &tau_f);
	fscanf(param, "uy_top %le\n", &uy_top);
	fscanf(param, "uy_bottom %le\n", &uy_bottom);
	fscanf(param, "yforce %le\n", &yforce);
	
	fscanf(param, "kappa %le\n", &kappa);
	fscanf(param, "rho %le\n", &rho);
	fscanf(param, "A_ldg %le\n", &A_ldg);
	fscanf(param, "U %le\n", &U);
	fscanf(param, "xi %le\n", &xi);
	fscanf(param, "Gamma_rot %le\n", &Gamma_rot);
	
	fscanf(param, "debug_on %d\n", &debug_on);
	
	fclose(param);
	
	itau_f = 1.0/tau_f;
	dt     = 1.0;
	qdt    = 0.5/(double)n_evol_Q;
	T      = third;	
	S      = getS();
}

//Print title information
void title()
{
	printf("-----------------------------------------\n");
	printf("--------3D LATTICE BOLTZMANN CODE--------\n");
	printf("--BASED ON ALGORITHM FROM YEOMANS ET AL--\n");
	printf("---AUTHORS: RUI ZHANG & TYLER ROBERTS----\n");
	printf("-----UNIVERSITY OF CHICAGO, JAN 2014-----\n");
	printf("-----------------------------------------\n");
}
