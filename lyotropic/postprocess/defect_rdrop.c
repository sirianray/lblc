/* 
 * To measure the space extension of defects in radial droplet
 * Inputs:  type_3d.out, phi_3d.out, Q_3d.out
 * Outputs: 
 *
 *
 * Author: Rui Zhang (Sirius)
 * Date:   Mar/1/2018
 * Copy Right Reserved
 * */

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "math.h"
#include <unistd.h>

int main(int argc, char *argv[]){
    FILE *param, *pfile, *tfile, *qfile, *ofile;
    int inverse=0, Nx, Ny, Nz, points, nlist, frame1=-2, frame2=-1, info;
    int iflag, ijunk, eof, frame, id, iarg, i, j, k, i2, j2, k2, id1, id2, id3, ilist, nx, ii;
	float a[9], w[3];
    double phic=0.5, Sc=0.3;
    double cx, cy, cz, cn, ir, exmax, exmin, eymax, eymin, ezmax, ezmin, x, y, z, w0, w1, w2, w3, phi0, rsq, wt, r, S0;
    char junk[256];

    int *type;
    double *phi, *S, *xlist, *ylist, *zlist;

    iarg=1;
    while (iarg<argc) {
        if (argv[iarg][0]=='c' || argv[iarg][0]=='C') {
            phic = atof(argv[iarg+1]);  // critical phi distinguishing two phases (0.5 by default)
            iarg+=2;
        } else if (argv[iarg][0]=='q' || argv[iarg][0]=='Q') {
            Sc = atof(argv[iarg+1]);    // critical S distinguishing defect
            iarg+=2;
        } else if (argv[iarg][0]=='s' || argv[iarg][0]=='S') {
            frame1 = atoi(argv[iarg+1]);// starting frame: frame1
            iarg+=2;
        } else if (argv[iarg][0]=='e' || argv[iarg][0]=='E') {
            frame2 = atoi(argv[iarg+1]);// ending frame: frame2
            iarg+=2;
        }
    }

    printf("phi_c = %f, S_c = %f\n",phic,Sc);

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
    type  = malloc(points*sizeof(int));
    phi   = malloc(points*sizeof(double));
    S     = malloc(points*sizeof(double));
    xlist = malloc(points*3*sizeof(double));
    ylist = malloc(points*3*sizeof(double));
    zlist = malloc(points*3*sizeof(double));
   
//  open input & output files
    tfile = fopen("type_3d.out","r");
    pfile = fopen("phi_3d.out","r");
    qfile = fopen("Q_3d.out","r");
    ofile = fopen("defect","w");

    eof   = 1;      // end of file flag
    frame = 1;
    while (eof>0 && (frame2<0 || frame<=frame2)) {
        for (id=0; id<points; id++) {
            if (tfile!=NULL) eof = fscanf(tfile, "%d\n", &type[id]);
            if (eof>0) eof = fscanf(pfile, "%lf\n", &phi[id]);
			if (eof>0) eof = fscanf(qfile,"%f %f %f %f %f\n",&a[0],&a[3],&a[6],&a[4],&a[7]);
            if (eof<=0) {
                printf("end of file\n");
                break;
            }

            if ((frame1<0 || frame>=frame1) && (frame2<0 || frame<=frame2) && phi[id]>phic) {
                a[8] = -1*(a[0] + a[4]);
                info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, w);

                if(info>0) {
	                printf("failure in eigenvalue routine\n");
                        exit(1);
                }

                S[id] = 0.5*3*w[2];
            }
        }

        if (eof>0 && (frame1<0 || frame>=frame1) && (frame2<0 || frame<=frame2)) {
            printf("frame %d\n",frame);

            nlist = 0;
            for (k=0; k<Nz; k++) {
                k2 = k+1;
                for (j=0; j<Ny; j++) {
                    j2 = j+1;
                    for (i=0; i<Nx; i++) {
                        i2 = i+1;
                        id = i + (j+k*Ny)*Nx;
                        S0 = S[id];
                        w0 = fabs(S0-Sc);
                        if (type[id]==-1 && phi[id]>phic) {
                            if (i2<Nx) {
                                id1 = id+1;
                                if (type[id1]==-1 && phi[id1]>phic && (S0-Sc)*(S[id1]-Sc)<=0) {
                                    w1 = fabs(S[id1]-Sc);
                                    if (w0+w1>0) {
                                        xlist[nlist] = (double)i + w0/(w0+w1);
                                    } else {
                                        xlist[nlist] = (double)i + 0.5;
                                    }
                                    ylist[nlist] = j;
                                    zlist[nlist] = k;
                                    nlist++;
                                }
                            }
                            if (j2<Ny) {
                                id2 = id + Nx;
                                if (type[id2]==-1 && phi[id2]>phic && (S0-Sc)*(S[id2]-Sc)<=0) {
                                    w2 = fabs(S[id2]-Sc);
                                    if (w0+w2>0) {
                                        ylist[nlist] = (double)j + w0/(w0+w2);
                                    } else {
                                        ylist[nlist] = (double)j + 0.5;
                                    }
                                    xlist[nlist] = i;
                                    zlist[nlist] = k;
                                    nlist++;
                                }
                            }
                            if (k2<Nz) {
                                id3 = id + Nx*Ny;
                                if (type[id3]==-1 && phi[id3]>phic && (S0-Sc)*(S[id3]-Sc)<=0) {
                                    w3 = fabs(phi[id3]-phic);
                                    if (w0+w3>0) {
                                        zlist[nlist] = (double)k + w0/(w0+w3);
                                    } else {
                                        zlist[nlist] = (double)k + 0.5;
                                    }
                                    xlist[nlist] = i;
                                    ylist[nlist] = j;
                                    nlist++;
                                }
                            }
                        }
                    }
                }
            }

            exmin = Nx;     // extension of tactoid
            exmax = 0;
            eymin = Ny;
            eymax = 0;
            ezmin = Nz;
            ezmax = 0;
            for (ilist=0; ilist<nlist; ilist++) {
                x = xlist[ilist];
                y = ylist[ilist];
                z = zlist[ilist];
                if (x>exmax) exmax=x;
                if (x<exmin) exmin=x;
                if (y>eymax) eymax=y;
                if (y<eymin) eymin=y;
                if (z>ezmax) ezmax=z;
                if (z<ezmin) ezmin=z;
            }

            if (nlist>0) {
                fprintf(ofile,"%f %f  %f %f  %f %f\n",exmin,exmax,eymin,eymax,ezmin,ezmax);
            } else {
                fprintf(ofile,"0. 0.  0. 0.  0. 0.\n");
            }
        }

        frame++;
    }

    if (tfile!=NULL) fclose(tfile);
    fclose(pfile);
    fclose(ofile);

    free(type);
    free(phi);
    free(S);
    free(xlist);
    free(ylist);
    free(zlist);
    return 0;
}
