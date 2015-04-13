#include "lb.h"

//Read in parameters from input FILE
void read_param()
{
	FILE *param;
	param = fopen("param.in", "r");

	fscanf(param, "newrun %d\n", &newrun);
	fscanf(param, "init_linear %d\n", &init_linear);
	fscanf(param, "simul_evol %d\n", &simul_evol);

	fscanf(param, "pbc_z %d\n", &pbc_z);
	fscanf(param, "wall_anchoring %d\n", &wall_anchoring);
	fscanf(param, "W_bot %lf\n", &W_bot);
	fscanf(param, "W_top %lf\n", &W_top);
	fscanf(param, "nbot %lf %lf %lf\n", &nbot[0], &nbot[1], &nbot[2]);
	fscanf(param, "ntop %lf %lf %lf\n", &ntop[0], &ntop[1], &ntop[2]);

	//Cavity dimentions
	fscanf(param, "cav_on %d\n", &cav_on); //Cavity switch (0=off, 1=on)
	fscanf(param, "cav_height %d\n", &cav_height); //A single cavity centered at Ly/2
	fscanf(param, "cav_width %d\n", &cav_width);
	fscanf(param, "ncavtop %lf %lf %lf\n", &ncavtop[0],  &ncavtop[1],  &ncavtop[2]  );
	fscanf(param, "ncavsides %lf %lf %lf\n",  &ncavsides[0],  &ncavsides[1],  &ncavsides[2]);

	fscanf(param, "rand_init %d\n", &rand_init);
	fscanf(param, "rand_seed %d\n", &rand_seed);
	fscanf(param, "ninit %lf %lf %lf\n", &ninit[0], &ninit[1], &ninit[2]);

	fscanf(param, "u_on %d\n", &u_on);
	fscanf(param, "Q_on %d\n", &Q_on);

	fscanf(param, "t_max %d\n", &t_max);
	fscanf(param, "t_v_on %d\n", &t_v_on);
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
	
	fscanf(param, "debug_mode %d\n", &debug_mode);
	fscanf(param, "movie %d\n", &movie);

	fclose(param);

	dt         = 1.0;
	itau_f     = 1.0/tau_f;
	itau_G     = 1.0/tau_G;
	T          = third;
	S = getS();

	if(cav_on==1 && uy_bottom!=0){
			  printf("ERROR: CANNOT HAVE BOTTOM WALL MOVING WITH CAVITY PRESENT!!!\n\n");
			  exit(1);
	}

			 
}


//Print title information
void title()
{
		  printf("-----------------------------------------\n");
		  printf("--3D LATTICE BOLTZMANN CODE--------------\n");
		  printf("--BASED ON ALGORITHM FROM YEOMANS ET AL.-\n");
		  printf("--AUTHORS: TYLER ROBERTS, RUI ZHANG------\n");
		  printf("--UNIVERSITY OF CHICAGO, JAN 2014--------\n");
		  printf("-----------------------------------------\n");
}
