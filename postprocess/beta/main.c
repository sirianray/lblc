/*
 * output Omega, beta for defect points
 * inputs:  param.in, Q_3d.out, type_3d.out
 * outputs: movie*.vtk
 *
 * Author: Rui Zhang (Sirius)
 * Date: Jun/14/2020
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

    int iarg, Nx=0, Ny=0, Nz=0, Np=0, S_flag=-1, ovtk=1, obeta=0, warning=1;
    int wall_x=0, wall_y=0, wall_z=0;
    int frame1=1, frame2=-1, eof, eof2, lines, id, i, j, k, l, itype;
    int in, jn, kn, ip, jp, kp, idp, im, jm, km;
    int ii, ijunk, iflag, info, jlist, ib;
    int ndpt, idpt, nnb=0;

    float a[9], w[3], b[9];
    double S0, Smax, Smin, Sth, rc=5., rc2=25., rd=3.1, rd2=9.61, zd=1.8, zd2=3.24;
    double dx, dy, dz, rsq, irsq, hLx, hLy, hLz, ir, dot, r2min;
    double tmpx, tmpy, tmpz, cph, sph, cth, sth, cth2, sth2, mx, my, mz;
    double M[6];
    complex tmpu, tmpv;

    char junk[256], fname[256], mfile[256];

    int *type, *nblist, *dpti, *dptc, *nbdid;
    float *nx, *ny, *nz, *S, *tx, *ty, *tz, *Ox, *Oy, *Oz, *Ux, *Uy, *Uz, *beta;
    float *dptx, *dpty, *dptz, *dptS, *dptM; 
    complex *dptA;

//  user-specified parameters
    iarg = 1;
    while ( iarg < argc ) {
        if (argv[iarg][0]=='r' && argv[iarg][1]=='c') {
            rc = atof(argv[iarg+1]);			    // rc for definition of neighbors solely for defects
            rc2= rc * rc;
            iarg+=2;
        } else if (argv[iarg][0]=='r' && argv[iarg][1]=='d') {
            rd = atof(argv[iarg+1]);			    // rd for defining bulk points neighboring to defects
            rd2= rd * rd;
            iarg+=2;
        } else if (argv[iarg][0]=='z') {
            zd = atof(argv[iarg+1]);			    // zd (half height of the cylinder) for calculating Omega
            zd2= zd * zd;
            iarg+=2;
        } else if (argv[iarg][0]=='s') {            // set threshold order parameter Sth
            Sth    = atof(argv[iarg+1]);
            S_flag = 1;
            iarg  += 2;
        } else if (argv[iarg][0]=='f' || argv[iarg][0]=='F') {
            frame1 = atoi(argv[iarg+1]);            // frames allowed to process
            iarg+=2;
            if (iarg<argc) {
                ijunk = atoi(argv[iarg]);
                if (ijunk>=frame1) {
                    frame2 = ijunk;
                    iarg++;
                }
            }
        } else if (argv[iarg][0]=='n' && argv[iarg][1]=='v') {
            ovtk = 0;                               // not ouput vtk files
            iarg++;
        } else if (argv[iarg][0]=='n' && argv[iarg][1]=='w') {
            warning = 0;                            // no warning
            iarg++;
        } else if (argv[iarg][0]=='b') {
            obeta = 1;                              // output beta files
            iarg++;
        } else iarg++;
    }

//  set parameters

//  find Nx, Ny, Nz from param.in or grid.out
    param = fopen("param.in","r");
    if (param != NULL) {
        iflag = 2; 
 //       if (S_flag == 0) iflag = 1;
        while ( iflag<3  && fgets (junk , 256 , param) != NULL ) {
            if (junk[0] == 'N' && junk[1] == 'x' && junk[2] == ' ') {
                sscanf(junk, "Nx Ny Nz %d %d %d\n", &Nx, &Ny, &Nz);	
                iflag++;
            }
            if (junk[0] == 'w' && junk[1] == 'a' && junk[2] == 'l') {
                sscanf(junk, "wall_x wall_y wall_z %d %d %d\n", &wall_x, &wall_y, &wall_z);	
                iflag++;
            }
//           if (S_flag==0 && junk[0] == 'U' && junk[1] == ' ') {
//               sscanf(junk, "U %lf\n", &S0);
//               if ( S0 >= 8./3. ) {               // equilibrium scalar order parameter
//                   S0 = 0.25 + 0.75 * sqrt(1. - 8./(3.*S0));
//                   S_flag = 1;
//               } else S_flag = -1;
//               iflag++;
//           }
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
    printf("box b.c. = [%d, %d, %d]\n", wall_x, wall_y, wall_z);
    printf("rc = %4.2f, rd = %4.2f, zd = %4.2f\n", rc, rd, zd);
    for (i=0; i<24; i++) printf("-");

//  half box sizes
    hLx   = 0.5 * (double)Nx;
    hLy   = 0.5 * (double)Ny;
    hLz   = 0.5 * (double)Nz;

//  allocation
    Np    = Nx * Ny * Nz;
    nx    = malloc(Np * sizeof(float));             // director
    ny    = malloc(Np * sizeof(float));
    nz    = malloc(Np * sizeof(float));
    S     = malloc(Np * sizeof(float));             // scalar order parameter
    type  = malloc(Np * sizeof(int));               // type (-1: bulk point)
    nblist= malloc(Np * sizeof(int));               // neighbor list (bulk points having defect neighbor)
    nbdid = malloc(Np * sizeof(int));               // nearest defect point id for neighbor list

    tx    = malloc(Np * sizeof(float));             // tangent of disclination curve (defect points)
    ty    = malloc(Np * sizeof(float));
    tz    = malloc(Np * sizeof(float));
    Ox    = malloc(Np * sizeof(float));             // Omega: rotation axis for director field (defect points)
    Oy    = malloc(Np * sizeof(float));
    Oz    = malloc(Np * sizeof(float));
    Ux    = malloc(Np * sizeof(float));             // U vector (orthogonal to t & Omega, defect points)
    Uy    = malloc(Np * sizeof(float));
    Uz    = malloc(Np * sizeof(float));
    beta  = malloc(Np * sizeof(float));             // beta angle between Omega & t (defect points & neighboring bulk points)

    dptx  = malloc(Np * sizeof(float));             // defect point position
    dpty  = malloc(Np * sizeof(float));
    dptz  = malloc(Np * sizeof(float));
    dptS  = malloc(Np * sizeof(float));             // defect point order parameter
    dpti  = malloc(Np * sizeof(int));               // defect point id
    dptc  = malloc(Np * sizeof(int));               // # of corresponding neighbor bulk points for defect points

    dptM  = malloc(Np * sizeof(float) * 6);         // auxiliary real matrix for defect points
    dptA  = malloc(Np * sizeof(complex) * 4);       // auxiliary complex matrix for defect points

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

        Smin = 100;
        Smax = -1.;

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

                if (Smin > S[id]) Smin = S[id];
                if (Smax < S[id]) Smax = S[id];
            }
            // post-read
			eof = fscanf(qfile,"%f %f %f %f %f\n",&a[0],&a[3],&a[6],&a[4],&a[7]);
            eof2= fscanf(tfile,"%d\n", &itype);
        }

        if (lines >= frame1) printf(":");
        printf("\n");
        fflush(stdout);

        if ( lines >= frame1) {
            if (S_flag < 0) Sth = 0.5 * (Smax + Smin);

            printf("\tSth = %5.3f\n",Sth);

            ndpt = 0;                                   // number of defect points

            if (Smin < Smax) {
                for (id = 0; id < Np; id++) {           // identify defect points
                    if (S[id] <= Sth && type[id]==-1) {
                        dptx[ndpt] = id%Nx;
                        dpty[ndpt] = (id/Nx)%Ny;
                        dptz[ndpt] = (id/Nx)/Ny;
                        dptS[ndpt] = S[id];
                        dpti[ndpt] = id;
                        ndpt++;
                    }
                }
            }

            printf("\tDefect point number = %d, %5.2f%%\n",ndpt,(double)ndpt/(double)Np*100.);
            if (ndpt*20 >= Np && warning != 0) {    // deal with too many defects
                printf("\tWarning: too many defect points. ");
                printf("Current S_threshold=%f, S_min=%f, S_mean=%f. Do you want to change S_threshold?\n", Sth, Smin, 0.5*(Smax+Smin));

                scanf("%s", junk);
                if (junk[0]=='Y' || junk[0]=='y') {
                    printf("\tplease enter new S_threshold");
                    scanf("%lf", &Sth);
                    printf("\n\tnow S_threshold=%f\n", Sth);

                    ndpt = 0;
                    for (id = 0; id < Np; id++) {
                        if (S[id] <= Sth && type[id]==-1) {
                            dpti[ndpt] = id;
                            dptx[ndpt] = id%Nx;
                            dpty[ndpt] = (id/Nx)%Ny;
                            dptz[ndpt] = (id/Nx)/Ny;
                            dptS[ndpt] = S[id];
                            ndpt++;
                        }
                    }
                    printf("\tAdjusted defect points = %d, %5.2f%%\n",ndpt, (double)ndpt/(double)Np*100.);
                }
            }

            for (id = 0; id < Np; id++) {           // initialization
                tx[id] = 0.;
                ty[id] = 0.;
                tz[id] = 0.;
                Ox[id] = 0.;
                Oy[id] = 0.;
                Oz[id] = 0.;
                Ux[id] = 0.;
                Uy[id] = 0.;
                Uz[id] = 0.;
                beta[id]=0.;
            }

            // initialize drdr matrix for defect points
            for (i = 0; i < ndpt * 6; i++) dptM[i] = 0.;

            for (i = 0; i < ndpt-1; i++) {          // cal dptM in order to cal tangent
                for (j = i+1; j < ndpt; j++) {
                    dx = dptx[i] - dptx[j];
                    dy = dpty[i] - dpty[j];
                    dz = dptz[i] - dptz[j];

                    if (wall_x != 0) {
                        if (dx > hLx) dx -= Nx;
                        if (dx <-hLx) dx += Nx;
                    }
                    if (wall_y != 0) {
                        if (dy > hLy) dy -= Ny;
                        if (dy <-hLy) dy += Ny;
                    }
                    if (wall_z != 0) {
                        if (dz > hLz) dz -= Nz;
                        if (dz <-hLz) dz += Nz;
                    }

                    rsq = dx * dx + dy * dy + dz * dz;
                    if (rsq <= rc2 && (rsq > zd2 || rc2 < zd2)) {
                        if (rsq > 0.) {
                            irsq = 1./rsq;
                            M[0] = dx * dx * irsq;
                            M[1] = dx * dy * irsq;
                            M[2] = dx * dz * irsq;
                            M[3] = dy * dy * irsq;
                            M[4] = dy * dz * irsq;
                            M[5] = dz * dz * irsq;

                            for (ii = 0; ii < 6; ii++) {
                                dptM[i * 6 + ii] += M[ii];
                                dptM[j * 6 + ii] += M[ii];
                            }
                        }
                    }
                }
            }


            for (i = 0; i < ndpt; i++) {            // cal tangent
                b[0] = dptM[i * 6];
                b[3] = dptM[i * 6 + 1];
                b[6] = dptM[i * 6 + 2];
                b[4] = dptM[i * 6 + 3];
                b[7] = dptM[i * 6 + 4];
                b[8] = dptM[i * 6 + 5];

			    info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, b, 3, w);
			    if (info > 0) {
			        printf("failure in eigenvalue routine\n");
			    	exit(1);
			    }

                id = dpti[i];
			    tx[id] = b[6];
			    ty[id] = b[7];
			    tz[id] = b[8];
            }

            // initialize nn matrix
            for (i = 0; i < ndpt * 6; i++) dptM[i] = 0.;

            nnb = 0;                                // number of bulk points having defect neighbors
            l   = rd + 1.;
            for (id = 0; id < Np; id++) {           // cal nblist (bulk points having defect neighbors)
                // sweep all non-defect bulk points
                if (S[id] > Sth && type[id] == -1) {
                    i    = id%Nx;
                    j    = (id/Nx)%Ny;
                    k    = (id/Nx)/Ny;
                    iflag= 0;                       // iflag = 1 if it has defect neighbors
                    r2min= rd * rd * 2.;            // distance to the closest defect neighbor
                    for (kn = k - l; kn <= k + l; kn++) {
                        kp = kn;
                        dz = kn - k;
                        if (kn < 0) {
                            if (wall_z == 0) kp += Nz;
                            else continue;
                        }
                        if (kn >= Nz) {
                            if (wall_z == 0) kp -= Nz;
                            else continue;
                        }
                        for (jn = j - l; jn <= j + l; jn++) {
                            jp = jn;
                            dy = jn - j;
                            if (jn < 0) {
                                if (wall_y == 0) jp += Ny;
                                else continue;
                            }
                            if (jn >= Ny) {
                                if (wall_y == 0) jp -= Ny;
                                else continue;
                            }
                            for (in = i - l; in <= i + l; in++) {
                                ip = in;
                                dx = in - i;
                                if (in < 0) {
                                    if (wall_x == 0) ip += Nx;
                                    else continue;
                                }
                                if (in >= Nx) {
                                    if (wall_x == 0) ip -= Nx;
                                    else continue;
                                }

                                idp = ip + (jp + kp * Ny) * Nx;
                                rsq = dx * dx + dy * dy + dz * dz;

                                if (S[idp] <= Sth && type[idp] == -1 && rsq <= rd2) {
                                    iflag = 1;
                                    if (rsq < r2min) {
                                        im = ip;    // nearest defect point
                                        jm = jp;
                                        km = kp;
                                    }
                                }
                            }
                        }
                    }

                    if (iflag == 1) {                // bulk point has a defect neighbor
                        nblist[nnb] = id;
                        nbdid[nnb]  = im + (jm + km * Ny) * Nx;
                        nnb++;
                        for (idpt = 0; idpt < ndpt; idpt++) {
                            idp = dpti[idpt];
                            ip  = idp%Nx;
                            jp  = (idp/Nx)%Ny;
                            kp  = (idp/Nx)/Ny;

                            dx  = im - ip;
                            dy  = jm - jp;
                            dz  = km - kp;
                            if (wall_x != 0) {
                                if (dx > hLx) dx -= Nx;
                                if (dx <-hLx) dx += Nx;
                            }
                            if (wall_y != 0) {
                                if (dy > hLy) dy -= Ny;
                                if (dy <-hLy) dy += Ny;
                            }
                            if (wall_z != 0) {
                                if (dz > hLz) dz -= Nz;
                                if (dz <-hLz) dz += Nz;
                            }
                            rsq = dx * dx + dy * dy + dz * dz;

                            dx  = i - ip;
                            dy  = j - jp;
                            dz  = k - kp;
                            if (wall_x != 0) {
                                if (dx > hLx) dx -= Nx;
                                if (dx <-hLx) dx += Nx;
                            }
                            if (wall_y != 0) {
                                if (dy > hLy) dy -= Ny;
                                if (dy <-hLy) dy += Ny;
                            }
                            if (wall_z != 0) {
                                if (dz > hLz) dz -= Nz;
                                if (dz <-hLz) dz += Nz;
                            }
                            dot = dx * tx[idp] + dy * ty[idp] + dz * tz[idp];

                            if (fabs(dot) <= zd && rsq <= zd2 * 4.) {
                                M[0] = nx[id] * nx[id];
                                M[1] = nx[id] * ny[id];
                                M[2] = nx[id] * nz[id];
                                M[3] = ny[id] * ny[id];
                                M[4] = ny[id] * nz[id];
                                M[5] = nz[id] * nz[id];

                                for (ii = 0; ii < 6; ii++) dptM[idpt * 6 + ii] += M[ii];
                            }
                        }
                    }
                }
            }

            for (idpt = 0; idpt < ndpt; idpt++) {   // cal Omega & U
                id   = dpti[idpt];
                b[0] = dptM[idpt * 6];
                b[3] = dptM[idpt * 6 + 1];
                b[6] = dptM[idpt * 6 + 2];
                b[4] = dptM[idpt * 6 + 3];
                b[7] = dptM[idpt * 6 + 4];
                b[8] = dptM[idpt * 6 + 5];

			    info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, b, 3, w);
			    if (info > 0) {
			        printf("failure in eigenvalue routine\n");
			    	exit(1);
			    }

                dot     = b[0] * tx[id] + b[1] * ty[id] + b[2] * tz[id];
                if (dot < 0.) for (ii = 0; ii < 3; ii++) b[ii] = -b[ii];
			    Ox[id]  = b[0];
			    Oy[id]  = b[1];
			    Oz[id]  = b[2];

                if (fabs(dot) < 1.) {               // find vector U orthonal to both tangent and Omega
                    tmpx= ty[id] * Oz[id] - tz[id] * Oy[id];
                    tmpy= tz[id] * Ox[id] - tx[id] * Oz[id];
                    tmpz= tx[id] * Oy[id] - ty[id] * Ox[id];
                } else if (fabs(tx[id]) < 1.) {
                    dot = tx[id];
                    tmpx= 1. - tx[id] * dot;
                    tmpy= 0. - ty[id] * dot;
                    tmpz= 0. - tz[id] * dot;
                } else {
                    dot = tx[id] + ty[id] + tz[id];
                    tmpx= 1. - tx[id] * dot;
                    tmpy= 1. - ty[id] * dot;
                    tmpz= 1. - tz[id] * dot;
                } 
                ir      = sqrt(1./(tmpx * tmpx + tmpy * tmpy + tmpz * tmpz));
                Ux[id]  = tmpx * ir;
                Uy[id]  = tmpy * ir;
                Uz[id]  = tmpz * ir;
            }

            // initialize auxiliary matrices and neighbor counter
            for (i = 0; i < ndpt * 4; i++) dptA[i] = 0.;
            for (i = 0; i < ndpt; i++) dptc[i] = 0;

            for (ib = 0; ib < nnb; ib++) {          // cal beta
                id = nblist[ib];                    // find bulk point id
                i  = id%Nx;
                j  = (id/Nx)%Ny;
                k  = (id/Nx)/Ny;
                im = nbdid[ib]%Nx;
                jm = (nbdid[ib]/Nx)%Ny;
                km = (nbdid[ib]/Nx)/Ny;
                for (idpt = 0; idpt < ndpt; idpt++) {
                    idp = dpti[idpt];
                    ip  = idp%Nx;
                    jp  = (idp/Nx)%Ny;
                    kp  = (idp/Nx)/Ny;

                    dx  = im - ip;
                    dy  = jm - jp;
                    dz  = km - kp;
                    if (wall_x != 0) {
                        if (dx > hLx) dx -= Nx;
                        if (dx <-hLx) dx += Nx;
                    }
                    if (wall_y != 0) {
                        if (dy > hLy) dy -= Ny;
                        if (dy <-hLy) dy += Ny;
                    }
                    if (wall_z != 0) {
                        if (dz > hLz) dz -= Nz;
                        if (dz <-hLz) dz += Nz;
                    }
                    rsq = dx * dx + dy * dy + dz * dz;

                    dx  = i - ip;
                    dy  = j - jp;
                    dz  = k - kp;
                    if (wall_x != 0) {
                        if (dx > hLx) dx -= Nx;
                        if (dx <-hLx) dx += Nx;
                    }
                    if (wall_y != 0) {
                        if (dy > hLy) dy -= Ny;
                        if (dy <-hLy) dy += Ny;
                    }
                    if (wall_z != 0) {
                        if (dz > hLz) dz -= Nz;
                        if (dz <-hLz) dz += Nz;
                    }
                    dot = dx * tx[idp] + dy * ty[idp] + dz * tz[idp];

                    if (fabs(dot) <= zd && rsq <= 4. * zd2) {
                        dx -= dot * tx[idp];        // make r = (dx, dy, dz) orthogonal to t
                        dy -= dot * ty[idp];
                        dz -= dot * tz[idp];
                        rsq = dx * dx + dy * dy + dz * dz;
                        if (rsq > 0.) {
                            ir  = sqrt(1. / rsq);
                            dx *= ir; 
                            dy *= ir;
                            dz *= ir;
                            dot = nx[id] * Ox[idp] + ny[id] * Oy[idp] + nz[id] * Oz[idp];
                            // convert n to m, orthogonal to Omega
                            mx  = nx[id] - dot * Ox[idp];
                            my  = ny[id] - dot * Oy[idp];
                            mz  = nz[id] - dot * Oz[idp];
                            rsq = mx * mx + my * my + mz * mz;
                            if (rsq > 0.) {
                                ir  = sqrt(1. / rsq);
                                mx *= ir;
                                my *= ir;
                                mz *= ir;

                                // phi: azimuthal angle
                                // cos(phi) = r . U
                                // sin(phi) = t . (U x r)
                                cph = dx * Ux[idp] + dy * Uy[idp] + dz * Uz[idp];
                                sph = tx[idp] * (Uy[idp] * dz - Uz[idp] * dy) - ty[idp] * (Ux[idp] * dz - Uz[idp] * dx) + tz[idp] * (Ux[idp] * dy - Uy[idp] * dx);
                                // th: angle corresponding to the director
                                // cos(th) = m . U
                                // sin(th) = O . (U x m)
                                cth = mx * Ux[idp] + my * Uy[idp] + mz * Uz[idp];
                                sth = Ox[idp] * (Uy[idp] * mz - Uz[idp] * my) - Oy[idp] * (Ux[idp] * mz - Uz[idp] * mx) + Oz[idp] * (Ux[idp] * my - Uy[idp] * mx);
                                cth2= 2. * cth * cth - 1.;
                                sth2= 2. * sth * cth;
                                // exp(I*(2*th - phi)) = cos(2th - phi) + I * sin(2th - phi)
                                tmpu= cth2 * cph + sth2 * sph + I * (sth2 * cph - cth2 * sph);
                                // exp(I*(2*th + phi)) = cos(2th + phi) + I * sin(2th + phi)
                                tmpv= cth2 * cph - sth2 * sph + I * (sth2 * cph + cth2 * sph);
                                dptA[idpt * 4]     += tmpu;
                                dptA[idpt * 4 + 1] += tmpu * conj(tmpu);
                                dptA[idpt * 4 + 2] += tmpv;
                                dptA[idpt * 4 + 3] += tmpv * conj(tmpv);
                                dptc[idpt]++;
                            }
                        }
                    }
                }
            }

            for (idpt = 0; idpt < ndpt; idpt++) {   // cal beta
                id   = dpti[idpt];
                tmpu = dptA[idpt * 4];
                tmpu = dptA[idpt * 4 + 1] * (double)dptc[idpt] - tmpu * conj(tmpu);
                tmpv = dptA[idpt * 4 + 2];
                tmpv = dptA[idpt * 4 + 3] * (double)dptc[idpt] - tmpv * conj(tmpv);

                dot  = tx[id] * Ox[id] + ty[id] * Oy[id] + tz[id] * Oz[id];
                if (dot > 1.) dot = 1.;
                if (cabs(tmpu) <= cabs(tmpv)) beta[id] = acos(dot);
                else beta[id] = PI - acos(dot);
            }

            for (ib = 0; ib < nnb; ib++) {          // neighbor bulk points assign nonzero beta
                id  = nblist[ib];
                beta[id] = beta[nbdid[ib]];
            }

            if (ovtk != 0) {                        // output vtk file
	            sprintf(mfile,"movie_%05d.vtk",lines);
                ofile = fopen(mfile,"w");

	            fprintf(ofile,"# vtk DataFile Version 3.0\n");
                fprintf(ofile,"vtk output\n");
                fprintf(ofile,"ASCII\n");
                fprintf(ofile,"DATASET STRUCTURED_GRID\n");

	            fprintf(ofile,"DIMENSIONS\t %d\t %d\t %d\t\n",Nx,Ny,Nz);
                fprintf(ofile,"POINTS\t %d\t float\n",Np);

                for (k = 0; k < Nz; k++) {
                    for (j = 0; j < Ny; j++) {
                        for (i = 0; i < Nx; i++) {
                            fprintf(ofile,"\t%0.2f\t%0.2f\t%0.2f\n", (double)i, (double)j, (double)k);
                        }
                    }
                }

                fprintf(ofile,"\n");

                fprintf(ofile,"POINT_DATA\t%d\n", Np);

                fprintf(ofile,"SCALARS type int 1\n");
	            fprintf(ofile,"LOOKUP_TABLE default\n");
                for(i = 0; i < Np; i++){
                    fprintf(ofile,"\t%d\n",type[i]);
                }
	            fprintf(ofile,"\n");

                fprintf(ofile,"SCALARS S_field float 1\n");
                fprintf(ofile,"LOOKUP_TABLE default\n");

	            for( i = 0; i < Np; i++) fprintf(ofile,"\t%lf\n",S[i]);

	            fprintf(ofile,"\n");
                fprintf(ofile,"SCALARS beta float 1\n");
                fprintf(ofile,"LOOKUP_TABLE default\n");

	            for( i = 0; i < Np; i++) fprintf(ofile,"\t%lf\n", beta[i]);

	            fprintf(ofile,"\n");
                fprintf(ofile,"SCALARS sin(beta) float 1\n");
                fprintf(ofile,"LOOKUP_TABLE default\n");

	            for( i = 0; i < Np; i++) fprintf(ofile,"\t%lf\n", sin(beta[i]));

	            fprintf(ofile,"\n");
                fprintf(ofile,"VECTORS directors float\n");

	            for(i = 0; i < Np; i++) fprintf(ofile,"\t%f\t%f\t%f\n", nx[i], ny[i], nz[i]);

	            fprintf(ofile,"\n");
                fprintf(ofile,"VECTORS tangents float\n");

	            for(i = 0; i < Np; i++) fprintf(ofile,"\t%f\t%f\t%f\n", tx[i], ty[i], tz[i]);

	            fprintf(ofile,"\n");
                fprintf(ofile,"VECTORS Omega float\n");

	            for(i = 0; i < Np; i++) fprintf(ofile,"\t%f\t%f\t%f\n", Ox[i], Oy[i], Oz[i]);

	            fprintf(ofile,"\n");
                fprintf(ofile,"VECTORS Uvec float\n");

	            for(i = 0; i < Np; i++) fprintf(ofile,"\t%f\t%f\t%f\n", Ux[i], Uy[i], Uz[i]);

	            printf("\tprinted vtk file %d\n",lines);
                fclose(ofile);
            }

            if (obeta != 0) {                       // output beta files
	            sprintf(mfile,"beta%05d",lines);
                ofile = fopen(mfile,"w");

                for (idpt = 0; idpt < ndpt; idpt++) {
                    id = dpti[idpt];
                    fprintf(ofile, "%f\n", beta[id]);
                }

	            printf("\tprinted beta file %d\n",lines);
                fclose(ofile);
            }
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
    free(nblist);
    free(nbdid);

    free(tx);
    free(ty);
    free(tz);
    free(Ox);
    free(Oy);
    free(Oz);
    free(Ux);
    free(Uy);
    free(Uz);
    free(beta);

    free(dptx);
    free(dpty);
    free(dptz);
    free(dptS);
    free(dpti);
    free(dptc);
    free(dptM);
    free(dptA);

    printf("done.\n");

    return 0;
}
