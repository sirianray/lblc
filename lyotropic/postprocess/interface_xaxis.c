/* 
 * To capture the shape of interface given x-axis is the symmetry axis
 * Inputs:  type_3d.out, phi_3d.out
 * Outputs: ishape
 *
 *
 * Author: Rui Zhang (Sirius)
 * Date:   Feb/28/2018
 * Copy Right Reserved
 * */

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "math.h"
#include <unistd.h>

int main(int argc, char *argv[]){
    FILE *param, *pfile, *tfile, *ofile;
    int inverse=0, Nx, Ny, Nz, points, nlist, Nd, frame1=-2, frame2=-1;
    int iflag, ijunk, eof, frame, id, iarg, i, j, k, i2, j2, k2, id1, id2, id3, ilist, nx, ii;
    double phic=0.5, del_x=1.;
    double cx, cy, cz, cn, ir, exmax, exmin, eymax, eymin, ezmax, ezmin, nxmax, nxmin, nymax, nymin, nzmax, nzmin, x, y, z, w0, w1, w2, w3, phi0, rsq, wt, r;
    double rmin, rmax, ax, ay, az, bx, by, bz;
    char junk[256];

    int *type, *shape_n;
    double *phi, *xlist, *ylist, *zlist, *shape_r;

    iarg=1;
    while (iarg<argc) {
        if (argv[iarg][0]=='i' || argv[iarg][0]=='I') {
            inverse = 1;                // inverse<=0: tactoid (phi>0.5 inside); >=1: inverse tactoid (phi<0.5 inside)
            iarg++;
        } else if (argv[iarg][0]=='c' || argv[iarg][0]=='C') {
            phic = atof(argv[iarg+1]);  // critical phi distinguishing two phases (0.5 by default)
            iarg+=2;
        } else if (argv[iarg][0]=='d' || argv[iarg][0]=='D') {
            del_x = atof(argv[iarg+1]); // del_x
            iarg+=2;
        } else if (argv[iarg][0]=='s' || argv[iarg][0]=='S') {
            frame1 = atoi(argv[iarg+1]);// starting frame: frame1
            iarg+=2;
        } else if (argv[iarg][0]=='e' || argv[iarg][0]=='E') {
            frame2 = atoi(argv[iarg+1]);// ending frame: frame2
            iarg+=2;
        }
    }

    printf("phi_c = %f, del_x=%f\n",phic,del_x);

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
    Nd     = (int)((double)Nx/del_x)+1;

//  allocation
    type  = malloc(points*sizeof(int));
    phi   = malloc(points*sizeof(double));
    xlist = malloc(points*3*sizeof(double));
    ylist = malloc(points*3*sizeof(double));
    zlist = malloc(points*3*sizeof(double));
    shape_n = malloc(Nd*sizeof(int));
    shape_r = malloc(Nd*sizeof(double));
   
//  open input & output files
    tfile = fopen("type_3d.out","r");
    pfile = fopen("phi_3d.out","r");
    ofile = fopen("ishape","w");

    eof   = 1;      // end of file flag
    frame = 1;
    while (eof>0 && (frame2<0 || frame<=frame2)) {
        cx = 0.;
        cy = 0.;
        cz = 0.;
        cn = 0.;
        nlist=0;
        for (id=0; id<points; id++) {
            if (tfile!=NULL) eof = fscanf(tfile, "%d\n", &type[id]);
            if (eof>0)  eof = fscanf(pfile, "%lf\n", &phi[id]);
            if (eof<=0) {
                printf("end of file\n");
                break;
            }

            if (type[id]==-1 && (phi[id]-phic)*((double)inverse-0.5)<=0) {
                i  = id%Nx;
                j  = (id/Nx)%Ny;
                k  = (id/Nx)/Ny;
                wt = fabs(phi[id]-phic);
                cx+= (double)i*wt;
                cy+= (double)j*wt;
                cz+= (double)k*wt;
                cn+= wt;
            }
        }

        if (eof>0 && (frame1<0 || frame>=frame1) && (frame2<0 || frame<=frame2)) {
            printf("frame %d\n",frame);

            if (cn>0.) {
                ir  = 1./cn;
                cx *= ir;
                cy *= ir;
                cz *= ir;
            }
            nlist = 0;
            for (k=0; k<Nz; k++) {
                k2 = k+1;
                for (j=0; j<Ny; j++) {
                    j2 = j+1;
                    for (i=0; i<Nx; i++) {
                        i2 = i+1;
                        id = i + (j+k*Ny)*Nx;
                        phi0 = phi[id];
                        w0 = fabs(phi0-phic);
                        if (type[id]==-1) {
                            if (i2<Nx) {
                                id1 = id+1;
                                if (type[id1]==-1 && (phi0-phic)*(phi[id1]-phic)<=0) {
                                    w1 = fabs(phi[id1]-phic);
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
                                if (type[id2]==-1 && (phi0-phic)*(phi[id2]-phic)<=0) {
                                    w2 = fabs(phi[id2]-phic);
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
                                if (type[id3]==-1 && (phi0-phic)*(phi[id3]-phic)<=0) {
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

            for (ii=0; ii<Nd; ii++) {
                shape_n[ii] = 0;
                shape_r[ii] = 0.;
            }
            for (ilist=0; ilist<nlist; ilist++) {
                x = xlist[ilist];
                y = ylist[ilist] - cy;
                z = zlist[ilist] - cz;
                nx= (int)(x/del_x);
                r = sqrt(y*y+z*z);
                shape_r[nx] += r;
                shape_n[nx] += 1;
            }

            for (ii=0; ii<Nd; ii++) if (shape_n[ii]>0) shape_r[ii]/=(double)shape_n[ii];

            if (nlist>0) {
                for (ii=0; ii<Nd; ii++) fprintf(ofile,"%f ",shape_r[ii]);
            } else {
                for (ii=0; ii<Nd; ii++) fprintf(ofile,"0. ");
            }
            fprintf(ofile,"\n");
        }

        frame++;
    }

    if (tfile!=NULL) fclose(tfile);
    fclose(pfile);
    fclose(ofile);

    free(type);
    free(phi);
    free(xlist);
    free(ylist);
    free(zlist);
    free(shape_n);
    free(shape_r);
    return 0;
}
