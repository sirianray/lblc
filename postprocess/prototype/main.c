/*
 * prototype postprocess program
 *
 * Author: Rui Zhang (Sirius)
 * Date: Oct/22/2020
 * Copyright reserved
 * */

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "math.h"
#include "complex.h"
#include <unistd.h>
#define PI 3.1415926535897931

int main(int argc, char *argv[]) {
    FILE *param, *tfile, *qfile, *ofile;

    int iarg, Nx=0, Ny=0, Nz=0, Np, S_flag;
    int frame1=1, frame2=-1, eof, eof2, lines, id, i, j, k, itype;
    int ijunk, iflag, info;

    float a[9], w[3];
    double S0;

    char junk[256], fname[256];

    int *type;
    float *nx, *ny, *nz, *S;

//  user-specified parameters
    iarg = 1;
    while ( iarg < argc ) {
        if (argv[iarg][0]=='f' || argv[iarg][0]=='F') {
            frame1 = atoi(argv[iarg+1]);			// frames allowed to process
            iarg+=2;
            if (iarg<argc) {
                ijunk = atoi(argv[iarg]);
                if (ijunk>=frame1) {
                    frame2 = ijunk;
                    iarg++;
                }
            }
        }
    }

//  set parameters

//  find Nx, Ny, Nz from param.in or grid.out; find U if initial S_flag=0
    param = fopen("param.in","r");
    if (param != NULL) {
        iflag = 1; 
        if (S_flag==0) iflag = 0;
        while ( iflag<2  && fgets (junk , 256 , param) != NULL ) {
            if (junk[0] == 'N' && junk[1] == 'x' && junk[2] == ' ') {
                sscanf(junk, "Nx Ny Nz %d %d %d\n", &Nx, &Ny, &Nz);	
                iflag++;
            }
            if (S_flag==0 && junk[0] == 'U' && junk[1] == ' ') {
                sscanf(junk, "U %lf\n", &S0);
                if ( S0 >= 8./3. ) {               // equilibrium scalar order parameter
                    S0 = 0.25 + 0.75 * sqrt(1. - 8./(3.*S0));
                    S_flag = 1;
                } else S_flag = -1;
                iflag++;
            }
        }
        fclose(param);
    } else {
        param = fopen("grid.out","r");
        if (param == NULL) {
            printf("error: system size info missing: check param.in or grid.out\n");
            exit(2);
        }
        fscanf(param, "Nx Ny Nz %d %d %d\n", &Nx, &Ny, &Nz);	
        fclose(param);
    }

//  display parameters
    printf("system parameters:\n");
    for (i=0; i<24; i++) printf("-");
    printf("\nbox size = [%d, %d, %d]\n", Nx, Ny, Nz);
    for (i=0; i<24; i++) printf("-");

//  allocation
    Np = Nx * Ny * Nz;
    nx = malloc(Np * sizeof(float));
    ny = malloc(Np * sizeof(float));
    nz = malloc(Np * sizeof(float));
    S  = malloc(Np * sizeof(float));
    type = malloc(Np * sizeof(int));

//  read Q_3d.out & type_3d.out, calculate image
    printf("\nreading:\n");
    qfile = fopen("Q_3d.out", "r");
    tfile = fopen("type_3d.out", "r");
    eof   = 0;      // end of file flag
    lines = 0;      // current frame

    // pre-read
	eof = fscanf(qfile,"%f %f %f %f %f\n",&a[0],&a[3],&a[6],&a[4],&a[7]);
    eof2= fscanf(tfile,"%d\n", &itype);

	while (eof>=0 && (lines<frame2 || frame2<0)) {
		lines++;
		printf("frame %d",lines);

//      reading Q & type
        for (id=0; id<Np; id++) {
            if (eof2> 0) type[id] = itype;
			if (eof > 0 && lines >= frame1 && type[id] == -1) {
				a[8] = -1*(a[0] + a[4]);
				info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, w);

				if(info>0) {
					printf("failure in eigenvalue routine\n");
					exit(1);
				}

				S[id]  = 1.5*w[2];
				nx[id] = a[6];
				ny[id] = a[7];
				nz[id] = a[8];
            }
            // post-read
			eof = fscanf(qfile,"%f %f %f %f %f\n",&a[0],&a[3],&a[6],&a[4],&a[7]);
            eof2= fscanf(tfile,"%d\n", &itype);
        }

        printf("\n");
        fflush(stdout);

        if ( lines >= frame1) {
        }
    }



//  close file
    fclose(qfile);
    fclose(tfile);

//  deallocation
    free(nx);
    free(ny);
    free(nz);
    free(S);
    free(type);

    printf("done.\n");

    return 0;
}
