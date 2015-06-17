#include "main.h"

void add_patch()
{
	FILE *param;
	int type_patch, i, j, id, jmin, jmax, ii;
	double plength, W_patch, x, y, z, q[6];

	if (patch_on!=0 && wall_on!=0 && Q_on!=0 && myid==0) {
		param = fopen("patch.in", "r");

		fscanf(param, "patch_length %lf\n", &plength);
		fscanf(param, "type W_patch %d %lf\n", &type_patch, &W_patch);
		fscanf(param, "n_patch %lf %lf %lf\n", &x, &y, &z);


		jmin = ceil(((double)Ny - 1.0 - plength) * 0.5);
		jmax =      ((double)Ny - 1.0 + plength) * 0.5;

		if (jmin<0)   jmin = 0;
		if (jmax>=Ny) jmax=Ny-1;

		for (j=jmin; j<=jmax; j++) {
			for (i=0; i<Nx; i++) {
				id = i + j*Nx;
				surf[id*10]  = type_patch;
				surf[id*10+1]= 2.*S_lc*S_lc*W_patch;
				ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
				for (ii=0; ii<5; ii++) surf[id*10+5+ii]=q[ii];
			}
		}

		fclose(param);
	}
}
