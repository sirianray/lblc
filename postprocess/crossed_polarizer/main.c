/*
 *
 * required inputs: 
 * default outputs:
 * optional outputs:
 *
 * Author: Rui Zhang (Sirius)
 * Date: Jun/11/2020
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

    int iarg, direc=0, Nx=0, Ny=0, Nz=0, Np, S_flag=-1;
    int frame1=1, frame2=-1, eof, eof2, lines, id, i, j, k, id1, id2, id3, N1, N2, N3, itype, idd;
    int ijunk, iflag, id1mask, id3mask, info;

    float a[9], w[3];

    double n_e=1.2, n_o=1., lambda=50., ilambda, vec_e[3]={0.}, S0=0.;
    double rsq, ir, ne, gi, sing, cosg, ai, sina, cosa, dot;
    double complex vec_p[2]={1.,0.}, vec_a[2]={0.,1.}, polx, poly, polxn, polyn;
    double complex e_ne, e_no, T[3]={0.};

    char junk[256], fname[256];

    int *type;
    float *nx, *ny, *nz, *S, *pol;

//  user-specified parameters
    iarg = 1;
    while ( iarg < argc ) {
        if (argv[iarg][0]=='x' || ((argv[iarg][0]=='+' || argv[iarg][0]=='-') && argv[iarg][1]=='x')) {
            direc = 1;                              // light source: +1 from x = +inf, -1 from x = -inf, etc.
            if (argv[iarg][0]=='-') direc = -1;
            iarg++;
        } else if (argv[iarg][0]=='y' || ((argv[iarg][0]=='+' || argv[iarg][0]=='-') && argv[iarg][1]=='y')) {
            direc = 2;
            if (argv[iarg][0]=='-') direc = -2;
            iarg++;
        } else if (argv[iarg][0]=='z' || ((argv[iarg][0]=='+' || argv[iarg][0]=='-') && argv[iarg][1]=='z')) {
            direc = 3;
            if (argv[iarg][0]=='-') direc = -3;
            iarg++;
        } else if (argv[iarg][0]=='l' || argv[iarg][0]=='L') {
            lambda = atof(argv[iarg+1]);            // incident light wavelength in lattice units (10 nm in most cases)
            iarg  += 2;
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
        } else if (argv[iarg][0]=='p') {            // polarizer (only real numbers)
            vec_p[0] = atof(argv[iarg+1]);
            vec_p[1] = atof(argv[iarg+2]);
            iarg += 3;
        } else if (argv[iarg][0]=='a') {            // analyzer (only real numbers)
            vec_a[0] = atof(argv[iarg+1]);
            vec_a[1] = atof(argv[iarg+2]);
            iarg += 3;
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

//  set parameters
    ilambda = 1./lambda;                            // 1/lambda
    if ( direc == 0) direc = 3;                     // default light source from z = +inf
    rsq = cabs(vec_p[0])*cabs(vec_p[0])+cabs(vec_p[1])*cabs(vec_p[1]);
    ir  = sqrt(1./rsq);
    vec_p[0] *= ir;
    vec_p[1] *= ir;
    rsq = cabs(vec_a[0])*cabs(vec_a[0])+cabs(vec_a[1])*cabs(vec_a[1]);
    ir  = sqrt(1./rsq);
    vec_a[0] *= ir;
    vec_a[1] *= ir;

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
    if (direc==3 || direc==-3) {
        N1 = Nx;
        N2 = Ny;
        N3 = Nz;
        vec_e[2] = direc/3;
    } else if (direc==2 || direc==-2) {
        N1 = Nx;
        N2 = Nz;
        N3 = Ny;
        vec_e[1] = direc/2;
    } else if (direc==1 || direc==-1) {
        N1 = Ny;
        N2 = Nz;
        N3 = Nx;
        vec_e[0] = direc;
    }

//  display parameters
    printf("system parameters:\n");
    for (i=0; i<24; i++) printf("-");
    printf("\nbox size = [%d, %d, %d]\n", Nx, Ny, Nz);
    printf("direc=[%2.0f, %2.0f, %2.0f], [N1, N2]=[%d, %d]\n", vec_e[0], vec_e[1], vec_e[2], N1, N2);
    printf("n_e=%f, n_o=%f\n", n_e, n_o);
    printf("incident light wavelength = %f\n", lambda);
    printf("polarizer=[%4.2f+%4.2fI, %4.2f+%4.2fI]\n", creal(vec_p[0]), cimag(vec_p[0]), creal(vec_p[1]), cimag(vec_p[1]));
    printf("analyzer =[%4.2f+%4.2fI, %4.2f+%4.2fI]\n", creal(vec_a[0]), cimag(vec_a[0]), creal(vec_a[1]), cimag(vec_a[1]));
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
    pol= malloc(N1 * N2 * sizeof(float));

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
//          diagram:
//          direc   id1     id2     id3
//          +3      i       j       k
//          -3      Nx-i-1  j       Nz-k-1
//          +2      Nx-i-1  k       j
//          -2      i       k       Ny-j-1
//          +1      j       k       i
//          -1      Ny-j-1  k       Nx-i-1
            for ( id1=0; id1<N1; id1++ ) {
                id1mask = id1;
                if (direc==-3 || direc==2 || direc==-1) id1mask = N1-1-id1;
                if (abs(direc)==1) j = id1mask;
                else i = id1mask;
                for ( id2=0; id2<N2; id2++ ) {
                    if (abs(direc)==3) j = id2;
                    else k = id2;

                    polx = vec_p[0];
                    poly = vec_p[1];
                    for ( id3=0; id3<N3; id3++ ) {
                        id3mask = id3;
                        if (direc<0) id3mask = N3-1-id3;
                        if (abs(direc)==3)      k = id3mask;
                        else if (abs(direc)==2) j = id3mask;
                        else if (abs(direc)==1) i = id3mask;

                        id = i + (j + k * Ny) * Nx;

                        if (type[id]==-1) {
                            dot = nx[id]*vec_e[0] + ny[id]*vec_e[1] + nz[id]*vec_e[2];
                            if (dot < 0.) {
                                nx[id] = -nx[id];
                                ny[id] = -ny[id];
                                nz[id] = -nz[id];
                                dot    = -dot;
                            }
                            if (dot >  1.) dot =  1.;

                            gi   = acos(dot);
                            sing = sin(gi);
                            cosg = cos(gi);

                            ne   = n_e;
                            if (S_flag>0) ne = n_o + (n_e-n_o) * S[id] / S0;
                            ne   = n_o*ne/sqrt(n_o*n_o*sing*sing + ne*ne*cosg*cosg);
                            e_ne = cos(ne*2.*PI*ilambda) + I*sin(ne*2.*PI*ilambda);
                            e_no = cos(n_o*2.*PI*ilambda)+ I*sin(n_o*2.*PI*ilambda);

                            switch (direc) {
                                case +3:
                                    ai = atan2(ny[id], nx[id]);
                                    break;
                                case -3:
                                    ai = atan2(ny[id],-nx[id]);
                                    break;
                                case +2:
                                    ai = atan2(nz[id],-nx[id]);
                                    break;
                                case -2:
                                    ai = atan2(nz[id], nx[id]);
                                    break;
                                case +1:
                                    ai = atan2(nz[id], ny[id]);
                                    break;
                                case -1:
                                    ai = atan2(nz[id],-ny[id]);
                                    break;
                                default:
                                    ai = atan2(ny[id], nx[id]);
                            }
                            sina = sin(ai);
                            cosa = cos(ai);

                            T[0] = e_ne*cosa*cosa + e_no*sina*sina;
                            T[1] = (e_ne - e_no) * sina * cosa;
                            T[2] = e_ne*sina*sina + e_no*cosa*cosa;

                            polxn= polx*T[0] + poly*T[1];
                            polyn= polx*T[1] + poly*T[2];
                            polx = polxn;
                            poly = polyn;
                        }
                    }
                    
                    idd          = id1 + id2 * N1;
                    pol[idd]     = cabs(vec_a[0] * polx + vec_a[1] * poly);
                    pol[idd]    *= pol[idd];
                }
            }

//          output polarization
            sprintf(fname,"pol%d",lines);
            ofile = fopen(fname, "w");
            for (i=0; i<N1; i++) {
                for (j=0; j<N2; j++) {
                    idd = i + j * N1;
                    if (j<N2-1) fprintf(ofile, "%f\t", pol[idd]);
                    else fprintf(ofile, "%f\n", pol[idd]);
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
    free(pol);

    printf("done.\n");

    return 0;
}
