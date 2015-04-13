#include "main.h"

//Read in parameters from input FILE
void read_param()
{
	double r2, ir;
	int i;	
	FILE *param;

	param = fopen("param.in", "r");
	
	fscanf(param, "newrun_on %d\n", &newrun_on);
	fscanf(param, "Nx Ny Nz %d %d %d\n", &Nx, &Ny, &Nz);	
	fscanf(param, "npar %d\n", &npar);
	fscanf(param, "wall_on %d\n", &wall_on);
	fscanf(param, "flow_on %d\n", &flow_on);
	fscanf(param, "Q_on %d\n", &Q_on);

	fscanf(param, "rand_init rand_seed %d %d\n", &rand_init, &rand_seed);
	
	fscanf(param, "t_max t_print t_write %d %d %d\n", &t_max, &t_print, &t_write);
	
	fscanf(param, "n_evol_Q %d\n", &n_evol_Q);
	
	fscanf(param, "type_bot W_bot n_bot %d %lf %lf %lf %lf\n", &type_bot, &W_bot, &n_bot[0], &n_bot[1], &n_bot[2]);
	fscanf(param, "type_top W_top n_top %d %lf %lf %lf %lf\n", &type_top, &W_top, &n_top[0], &n_top[1], &n_top[2]);
	
	fscanf(param, "tau_f %le\n", &tau_f);
	fscanf(param, "uy_top uy_bottom %le %le\n", &uy_top, &uy_bottom);
	fscanf(param, "yforce %le\n", &yforce);
	
	fscanf(param, "kappa %le\n", &kappa);
	fscanf(param, "rho %le\n", &rho);
	fscanf(param, "A_ldg %le\n", &A_ldg);
	fscanf(param, "U %le\n", &U_lc);
	fscanf(param, "xi %le\n", &xi);
	fscanf(param, "Gamma_rot %le\n", &Gamma_rot);
	
	fscanf(param, "debug_on %d\n",&debug_on);
	
	fclose(param);

	points = Nx*Ny*Nz;
	if (wall_on!=0) {
		bulk0  = Nx*Ny;
		qpoints= Nx*Ny*(Nz+2);
		nsurf  = Nx*Ny*2;
	} else {
		bulk0  = 0;
		qpoints= points;
		nsurf  = 0;
	}
	
	itau_f = 1.0/tau_f;
	qdt    = 1/(double)n_evol_Q;
	S_lc   = getS(U_lc);	
	xi1    = 0.5*(xi+1.0);
	xi2    = 0.5*(xi-1.0);
	
	if (flow_on!=0) uconverge=0;
	if (Q_on!=0)    qconverge=0;

        r2 = n_top[0]*n_top[0]+n_top[1]*n_top[1]+n_top[2]*n_top[2];
        if (r2>0 && r2!=1) {
                ir = 1.0/sqrt(r2);
                for (i=0; i<3; i++) {
                        n_top[i] *= ir;
                }
        } else if (r2==0) {
                printf("Top wall preferred orientation wrong\n");
                exit(-1);
        }

        r2 = n_bot[0]*n_bot[0]+n_bot[1]*n_bot[1]+n_bot[2]*n_bot[2];
        if (r2>0 && r2!=1) {  
                ir = 1.0/sqrt(r2);   
                for (i=0; i<3; i++) {  
                        n_bot[i] *= ir;
                }
        } else if (r2==0) {  
                printf("Bottom wall preferred orientation wrong\n");
                exit(-1);
        }
}

//Print title information
void title()
{
	printf("-----------------------------------------\n");
	printf("------3D LATTICE BOLTZMANN MPI CODE------\n");
	printf("--BASED ON ALGORITHM FROM YEOMANS ET AL--\n");
	printf("---AUTHORS: RUI ZHANG & TYLER ROBERTS----\n");
	printf("-----UNIVERSITY OF CHICAGO, JAN 2014-----\n");
	printf("-----------------------------------------\n");
}
