#include "lb.h"

void alloc(){

		  int i, points = Nx*Ny*Nz;

		  Rho      = malloc(points*sizeof(double));
		  cav_flag = malloc(points*sizeof(double));
		  for(i=0;i<points;i++) {
					 Rho[i]  = 0;
					 cav_flag[i] = 0;
		  }

		  //[0->4 = 00, 01, 02, 11, 12]
		  Q1    = malloc(5*points*sizeof(double));
		  Q2    = malloc(5*points*sizeof(double));
		  convQ = malloc(5*points*sizeof(double));
		  H     = malloc(5*points*sizeof(double));
		  for(i=0;i<(5*points);i++){
					 Q1[i] = 0;
					 Q2[i] = 0;
					 H[i]  = 0;
					 convQ[i] = 0;
		  }

		  //[0->2 = ux, uy, uz]
		  u    = malloc(3*points*sizeof(double));
		  dtau = malloc(3*points*sizeof(double));
		  for(i=0;i<(3*points);i++) {
					 u[i] = 0;
					 dtau[i] = 0;
		  }

		  f1   = malloc(15*points*sizeof(double));
		  f2   = malloc(15*points*sizeof(double));
		  f_eq = malloc(15*points*sizeof(double));
		  Cf   = malloc(15*points*sizeof(double));
		  p    = malloc(15*points*sizeof(double));
		  for(i=0;i<(15*points);i++){
					 f1[i] = 0;
					 f2[i] = 0;
					 f_eq[i] = 0;
					 Cf[i] = 0;
					 p[i] = 0;
		  }

		  //In order [00 01 02; 10 11 12, 20 21 22]
		  W       = malloc(9*points*sizeof(double));
		  tau     = malloc(9*points*sizeof(double));
		  for(i=0;i<(9*points);i++) {
					 W[i] = 0;
					 tau[i] = 0;
		  }

		  sigma   = malloc(6*points*sizeof(double));
		  sigma_q   = malloc(6*points*sizeof(double));
		  for(i=0;i<(6*points);i++) {
					 sigma[i] = 0;
					 sigma_q[i] = 0;
		  }

		  //0->bottom, 1-> top. 0->5 00, 01, 02, 11, 12
		  Qsurfbot   = malloc(5*Nx*Ny*sizeof(double));
		  Qsurftop   = malloc(5*Nx*Ny*sizeof(double));
		  dQdxnu_top = malloc(5*Nx*Ny*sizeof(double));
		  dQdxnu_bot = malloc(5*Nx*Ny*sizeof(double));
		  for(i=0;i<(5*Nx*Ny);i++) {
					 Qsurfbot[i] = 0;
					 Qsurftop[i] = 0;
					 dQdxnu_top[i] = 0;
					 dQdxnu_bot[i] = 0;
		  }

		  printf("Allocated variables\n");

}

void free_all(){

		  free(Rho      );
		  free(cav_flag );

		  free(Q1 );
		  free(Q2 );
		  free(convQ );

		  free(u );
		  free(dtau );

		  free(f1   );
		  free(f2   );
		  free(f_eq );
		  free(Cf   );
		  free(p    );

		  free(H       );
		  free(W       );
		  free(sigma   );
		  free(sigma_q );
		  free(tau     );

		  free(Qsurftop  );
		  free(Qsurfbot  );
		  free(dQdxnu_top );
		  free(dQdxnu_bot );
}

