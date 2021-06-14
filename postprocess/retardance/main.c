/*
 *
 * required inputs: 
 * default outputs:
 * optional outputs:
 *
 * Author: Rui Zhang (Sirius)
 * Date: Jun/11/2020 (crossed_polarizer)
 * Date: Feb/27/2021
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

    int iarg, direc=1, Nx=0, Ny=0, Nz=0, Np, S_flag=-1;
    int frame1=1, frame2=-1, eof, eof2, lines, id, i, j, k, id1, id2, id3, N1, N2, N3, itype, idd;
    int ijunk, iflag, id1mask, id2mask, id3mask, info;

    float a[9], w[3];

    double n_e=1.2, n_o=1., lambda, ilambda, vec_e[3]={0.}, S0=0.;
    double rsq, ir, ne, gi, sing2, cosg2, ai, sina, cosa, dot, dis1, dis2;
    double complex m0, m1, m2, m3, ctmp, lmb1, lmb2;
    double complex e_ne, e_no, T[3]={0.}, M[2][2];

    char junk[256], fname[256];

    int *type;
    float *nx, *ny, *nz, *S;
    double *rtd;

//  user-specified parameters
    iarg = 1;
    while ( iarg < argc ) {
        if (argv[iarg][0]=='x') {
            direc = 1;                              // light source, in x by default
            iarg++;
        } else if (argv[iarg][0]=='y') {
            direc = 2;
            iarg++;
        } else if (argv[iarg][0]=='z') {
            direc = 3;
            iarg++;
        } else if (argv[iarg][0]=='s') {            // order parameter is now considered (by default, S_flag<0)
            S_flag = 0;                             // will read order parameter (essentially U) from param.in
            iarg++;
        } else if (argv[iarg][0]=='S') {
            S_flag = 1;                             // S0 is assigned by user
            S0 = atof(argv[iarg+1]);
            iarg  += 2;
        } else if (argv[iarg][0]=='n' && argv[iarg][1]=='o') {
            n_o  = atof(argv[iarg+1]);              // refraction index for ordinary direction
            iarg+= 2;
        } else if (argv[iarg][0]=='n' && argv[iarg][1]=='e') {
            n_e  = atof(argv[iarg+1]);              // refraction index for extraordinary direction
            iarg+= 2;
        } else if (argv[iarg][0]=='f' || argv[iarg][0]=='F') {
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

    printf("reading param.in...\n");

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

//  [N1, N2] is the projection plane (coordinate system for crossed polarization image)
//  N3 is the depth the light travels
    if (direc==3) {
        N1 = Nx;
        N2 = Ny;
        N3 = Nz;
        vec_e[2] = 1.;
        lambda = 4 * Nz;
    } else if (direc==2) {
        N1 = Nx;
        N2 = Nz;
        N3 = Ny;
        vec_e[1] = 1.;
        lambda = 4 * Ny;
    } else if (direc==1) {
        N1 = Ny;
        N2 = Nz;
        N3 = Nx;
        vec_e[0] = 1.;
        lambda = 4 * Nx;
    }
    ilambda = 1./(double)lambda;

//  display parameters
    printf("system parameters:\n");
    for (i=0; i<24; i++) printf("-");
    printf("\nbox size = [%d, %d, %d]\n", Nx, Ny, Nz);
    printf("direc=[%2.0f, %2.0f, %2.0f], [N1, N2]=[%d, %d]\n", vec_e[0], vec_e[1], vec_e[2], N1, N2);
    printf("n_e=%f, n_o=%f\n", n_e, n_o);
    if ( S0 > 0.) printf("order parameter S0 = %f\n", S0);
    else printf("order parameter not considered (by default)\n");
    for (i=0; i<24; i++) printf("-");

//  allocation
    Np = Nx * Ny * Nz;
    nx = malloc(Np * sizeof(float));
    ny = malloc(Np * sizeof(float));
    nz = malloc(Np * sizeof(float));
    S  = malloc(Np * sizeof(float));
    type = malloc(Np * sizeof(int));
    rtd= malloc(N1 * N2 * sizeof(double));

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
            for ( id1=0; id1<N1; id1++ ) {
                id1mask = id1;
                if (direc == 1) j = id1mask;
                else i = id1mask;
                for ( id2=0; id2<N2; id2++ ) {
                    id2mask = id2;
                    if (direc == 3) j = id2mask;
                    else k = id2mask;

                    M[0][0] = 1.;
                    M[0][1] = 0.;
                    M[1][0] = 0.;
                    M[1][1] = 1.;
                    for ( id3=0; id3<N3; id3++ ) {
                        id3mask = id3;
                        if (direc == 1) i = id3mask;
                        else if (direc == 2) j = id3mask;
                        else if (direc == 3) k = id3mask;

                        id = i + (j + k * Ny) * Nx;

                        if (type[id]==-1) {
                            dot = nx[id]*vec_e[0] + ny[id]*vec_e[1] + nz[id]*vec_e[2];
                            if (dot < 0.) {
                                nx[id] = -nx[id];
                                ny[id] = -ny[id];
                                nz[id] = -nz[id];
                                dot    = -dot;
                            }
                            cosg2 = dot*dot;
                            sing2 = 1. - cosg2;

                            ne   = n_e;
                            if (S_flag>0) ne = n_o + (n_e-n_o) * S[id] / S0;
                            ne   = n_o*ne/sqrt(n_o*n_o*sing2 + ne*ne*cosg2);
                            e_ne = cos(ne* 2.*PI*ilambda) + I*sin(ne *2.*PI*ilambda);
                            e_no = cos(n_o*2.*PI*ilambda) + I*sin(n_o*2.*PI*ilambda);

                            switch (direc) {
                                case +3:
                                    ai = atan2(ny[id], nx[id]);
                                    break;
                                case +2:
                                    ai = atan2(nz[id],-nx[id]);
                                    break;
                                case +1:
                                    ai = atan2(nz[id], ny[id]);
                                    break;
                                default:
                                    ai = atan2(ny[id], nx[id]);
                            }
                            sina = sin(ai);
                            cosa = cos(ai);

                            T[0] = e_ne*cosa*cosa + e_no*sina*sina;
                            T[1] = (e_ne - e_no) * sina * cosa;
                            T[2] = e_ne*sina*sina + e_no*cosa*cosa;

                            m0   = T[0]*M[0][0] + T[1]*M[1][0];
                            m1   = T[0]*M[0][1] + T[1]*M[1][1];
                            m2   = T[1]*M[0][0] + T[2]*M[1][0];
                            m3   = T[1]*M[0][1] + T[2]*M[1][1];
                            M[0][0] = m0;
                            M[0][1] = m1;
                            M[1][0] = m2;
                            M[1][1] = m3;
                        }
                    }

                    ctmp = (M[0][0] - M[1][1])*(M[0][0] - M[1][1]) + 4.*M[0][1]*M[1][0];
                    ctmp = csqrt(ctmp);
                    lmb1 = 0.5 * ( M[0][0] + M[1][1] + ctmp);
                    dis1 = carg(lmb1);
                    lmb2 = 0.5 * ( M[0][0] + M[1][1] - ctmp);
                    dis2 = carg(lmb2);

                    
                    idd  = id1 + id2 * N1;
                    rtd[idd] = fabs(dis1 - dis2)*lambda/(2.*PI); 
                }
            }

//          output retardance
            sprintf(fname,"retard%d",lines);
            ofile = fopen(fname, "w");
            for (i=0; i<N1; i++) {
                for (j=0; j<N2; j++) {
                    idd = i + j * N1;
                    if (j<N2-1) fprintf(ofile, "%f\t", rtd[idd]);
                    else fprintf(ofile, "%f\n", rtd[idd]);
                }
            }
            fclose(ofile);
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
    free(rtd);

    printf("done.\n");

    return 0;
}
