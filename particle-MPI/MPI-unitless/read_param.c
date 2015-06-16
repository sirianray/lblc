#include "main.h"

//Read in parameters from input FILE
void read_param()
{
	int ipar, i;
	double x, y, z, r, Wt, W, r2, ir;
	FILE *param;

	param = fopen("param.in", "r");
	
	fscanf(param, "newrun_on %d\n", &newrun_on);
	fscanf(param, "Nx Ny Nz %d %d %d\n", &Nx, &Ny, &Nz);	
	fscanf(param, "npar %d\n", &npar);
	fscanf(param, "patch_on %d\n", &patch_on);
	fscanf(param, "wall_on %d\n", &wall_on);
	fscanf(param, "flow_on %d\n", &flow_on);
	fscanf(param, "Q_on %d\n", &Q_on);

	fscanf(param, "rand_init rand_seed %d %d\n", &rand_init, &rand_seed);
	
	fscanf(param, "t_max t_print t_write %d %d %d\n", &t_max, &t_print, &t_write);
	
	fscanf(param, "n_evol_Q %d\n", &n_evol_Q);
	
	fscanf(param, "type_bot W_bot n_bot %d %lf %lf %lf %lf\n", &type_bot, &W_bot, &n_bot[0], &n_bot[1], &n_bot[2]);
	fscanf(param, "type_top W_top n_top %d %lf %lf %lf %lf\n", &type_top, &W_top, &n_top[0], &n_top[1], &n_top[2]);
	
	fscanf(param, "tau_f %le\n", &tau_f);
        fscanf(param, "ux_bot uy_bot %le %le\n", &ux_bot, &uy_bot);
	fscanf(param, "ux_top uy_top %le %le\n", &ux_top, &uy_top);
	fscanf(param, "xforce yforce %le %le\n", &xforce, &yforce);
	
	fscanf(param, "l0 %le\n", &l0);
	fscanf(param, "tau %le\n", &tau);
	fscanf(param, "nu %le\n", &nu);
	fscanf(param, "rho %le\n", &rho);
	fscanf(param, "U %le\n", &U_lc);
	fscanf(param, "xi %le\n", &xi);
	
	fscanf(param, "debug_on %d\n",&debug_on);
	
	fclose(param);

	if (npar>0) param = fopen("ppos.in","r");

	points = Nx*Ny*Nz;
	point=ceil((double)points/(double)numprocs);
	lpoint=point;
	if(myid==numprocs-1)lpoint=points-(numprocs-1)*point;

	if (wall_on!=0) {
		bulk0  = Nx*Ny;
		qpoints= Nx*Ny*(Nz+2);
		nsurf  = Nx*Ny*2;
		nodes  = Nx*Ny*2;
	} else {
		bulk0  = 0;
		qpoints= points;
		nsurf  = 0;
		nodes  = 0;
	}

	for (ipar=0; ipar<npar; ipar++) {
		fscanf(param, "%le %le %le %le %le %le\n", &x, &y, &z, &r, &Wt, &W);
		if (r<4.0) {
			nodes += ceil(r*r*31.0);
		} else if (r<8.0) {
                	nodes += ceil(r*r*22.0);
		} else if (r<14.0) {
                	nodes += ceil(r*r*20.0);
		} else {
                        nodes += ceil(r*r*19.5);
		}
	}
	if (npar>0) fclose(param);

	node  =ceil((double)nodes/(double)numprocs);
	qpoint=ceil((double)qpoints/(double)numprocs);

	if(myid==0) printf("nodes, node=%d %d\n",nodes,node);

	itau_f = 1.0/tau_f;
	qdt    = 1/(double)n_evol_Q;
	S_lc   = getS(U_lc);	
	Fld0   = getF0(U_lc);
	xi1    = 0.5*(xi+1.0);
	xi2    = 0.5*(xi-1.0);

	strcf  = nu *tau*tau/(l0*l0);

	Gamma_rot = 2.*S_lc*S_lc*tau;
	kappa  = 1./(l0*l0); 
	
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
	MPI_Barrier(MPI_COMM_WORLD);	
}

//Print title information
void title()
{
	if (myid == 0){
		printf("-----------------------------------------\n");
		printf("------3D LATTICE BOLTZMANN MPI CODE------\n");
		printf("--BASED ON ALGORITHM FROM YEOMANS ET AL--\n");
		printf("---AUTHORS: RUI ZHANG & TYLER ROBERTS----\n");
		printf("-----UNIVERSITY OF CHICAGO, JAN 2014-----\n");
		printf("-----------------------------------------\n");
	}
}
