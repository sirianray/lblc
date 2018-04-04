/* 
 * To measure vmax & sum_vsq in the system
 * Inputs:  type_3d.out, u_3d.out
 * Outputs: vgrid 
 *
 *
 * Author: Rui Zhang (Sirius)
 * Date:   Mar/1/2018
 * Copy Right Reserved
 * */

#include <stdio.h>
#include <stdlib.h>
//#include "mkl.h"
#include "math.h"
#include <unistd.h>

int main(int argc, char *argv[]){
    FILE *param, *ufile, *ofile;
    int inverse=0, Nx, Ny, Nz, points, nlist, frame1=-2, frame2=-1, info;
    int iflag, ijunk, eof, frame, id, iarg, i, j, k, i2, j2, k2, id1, id2, id3, ilist, iu;
	float a[9], w[3];
    double phic=0.5, vmsq, vmax;
    double cx, cy, cz, cn, ir, exmax, exmin, eymax, eymin, ezmax, ezmin, x, y, z, w0, w1, w2, w3, phi0, rsq, wt, r, S0, usq;
    char junk[256];

    int *type;
    double *u;

    iarg=1;
    while (iarg<argc) {
        if (argv[iarg][0]=='s' || argv[iarg][0]=='S') {
            frame1 = atoi(argv[iarg+1]);// starting frame: frame1
            iarg+=2;
        } else if (argv[iarg][0]=='e' || argv[iarg][0]=='E') {
            frame2 = atoi(argv[iarg+1]);// ending frame: frame2
            iarg+=2;
        }
    }

//  capture boundary condition from param.in
    param = fopen("param.in","r");
	iflag = 1;		// skip info lines which start with #
	while ( iflag==1 && fgets (junk , 256 , param) != NULL ) {
		if(junk[0]=='#') {
			iflag=1;
		} else {
			iflag=0;
		}
	}
	sscanf(junk,  "newrun_on %d\n", &ijunk);
	fscanf(param, "Nx Ny Nz %d %d %d\n", &Nx, &Ny, &Nz);	
//	fscanf(param, "npar %d\n", &npar);
//	fscanf(param, "patch_on %d\n", &ijunk);
//	fscanf(param, "wall_x wall_y wall_z %d %d %d\n", &wall_x, &wall_y, &wall_z);
    fclose(param);
    points = Nx*Ny*Nz;

//  allocation
    u   = malloc(points*3*sizeof(double));
   
//  open input & output files
    ufile = fopen("u_3d.out","r");
    ofile = fopen("vgrid","w");

    eof   = 1;      // end of file flag
    frame = 1;
    while (eof>0 && (frame2<0 || frame<=frame2)) {
        for (id=0; id<points; id++) {
            if (eof>0) eof = fscanf(ufile, "%lf %lf %lf\n",&u[id*3],&u[id*3+1],&u[id*3+2]);
            if (eof<=0) {
                printf("end of file\n");
                break;
            }
        }

        if (eof>0 && (frame1<0 || frame>=frame1) && (frame2<0 || frame<=frame2)) {
            printf("frame %d\n",frame);

            vmsq = 0.;
            vmax = 0.;
            for (id=0; id<points; id++) {
                iu = id*3;
                usq   = u[iu]*u[iu] + u[iu+1]*u[iu+1] + u[iu+2]*u[iu+2];
                vmsq += usq;
                if (usq>vmax) vmax = usq;
            }
            vmsq = sqrt(vmsq);
            vmax = sqrt(vmax);

            fprintf(ofile,"%e %e\n",vmax,vmsq);
        }

        frame++;
    }

    fclose(ufile);
    fclose(ofile);

    free(u);
    return 0;
}
