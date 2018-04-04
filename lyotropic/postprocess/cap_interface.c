/* 
 * To measure the elongation, neck radius of coalescing droplets
 * Inputs:  type_3d.out, phi_3d.out
 * Outputs: coalesce
 *
 * elements in coalesce: cx cy cz  ex_x ex_y ex_z  neck_x neck_y neck_z  rmax ax ay az  rmin bx by bz 
 *
 * c?: droplet center of mass (COM)
 * ex_?:   extension of interface in 3 directions
 * neck_?: neck diameter in 3 directions
 * rmax:  maximum distance of interfacial points to COM (semi-long-axis)
 * a?:    vector associated with rmax point
 * rmin:  minimum distance of interfacial points to COM (probably neck radius)
 * b?:    vector associated with rmin point
 *
 * Author: Rui Zhang (Sirius)
 * Date:   Feb/25/2018
 * Copy Right Reserved
 * */

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "math.h"
#include <unistd.h>

int main(int argc, char *argv[]){
    FILE *param, *pfile, *tfile, *ofile;
    int inverse=0, Nx, Ny, Nz, points, nlist;
    int iflag, ijunk, eof, frame, id, iarg, i, j, k, i2, j2, k2, id1, id2, id3, ilist;
    double phic=0.5;
    double cx, cy, cz, cn, ir, exmax, exmin, eymax, eymin, ezmax, ezmin, nxmax, nxmin, nymax, nymin, nzmax, nzmin, x, y, z, w0, w1, w2, w3, phi0, rsq, wt;
    double rmin, rmax, ax, ay, az, bx, by, bz;
    char junk[256];

    int *type;
    double *phi, *xlist, *ylist, *zlist;

    iarg=1;
    while (iarg<argc) {
        if (argv[iarg][0]=='i' || argv[iarg][0]=='I') {
            inverse = 1;                // inverse<=0: tactoid (phi>0.5 inside); >=1: inverse tactoid (phi<0.5 inside)
            iarg++;
        } else if (argv[iarg][0]=='c' || argv[iarg][0]=='C') {
            phic = atof(argv[iarg+1]);  // critical phi distinguishing two phases (0.5 by default)
            iarg+=2;
        }
    }

    printf("phi_c = %f\n",phic);

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
    xlist = malloc(points*3*sizeof(double));
    ylist = malloc(points*3*sizeof(double));
    zlist = malloc(points*3*sizeof(double));
   
//  open input & output files
    tfile = fopen("type_3d.out","r");
    pfile = fopen("phi_3d.out","r");
    ofile = fopen("coalesce","w");

    eof   = 1;      // end of file flag
    frame = 1;
    while (eof>0) {
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

        if (eof>0) {
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
                                    //if (frame==79) printf("%lf %lf %lf\n",xlist[nlist],ylist[nlist],zlist[nlist]);
                                    //if (frame==19) printf("%d %d %d  %f %f\n",i,j,k,phi0,phi[id1]);
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
                                    //if (frame==79) printf("%lf %lf %lf\n",xlist[nlist],ylist[nlist],zlist[nlist]);
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
                                    //if (frame==79) printf("%lf %lf %lf\n",xlist[nlist],ylist[nlist],zlist[nlist]);
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
            nxmin = Nx;     // neck of tactoid
            nxmax = 0;
            nymin = Ny;
            nymax = 0;
            nzmin = Nz;
            nzmax = 0;
            rmin  = Nx;
            rmax  = 0.;
            ax    = 0.;
            ay    = 0.;
            az    = 0.;
            bx    = 0.;
            by    = 0.;
            bz    = 0.;
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
                if ((y-cy)*(y-cy)+(z-cz)*(z-cz)<=1.) {
                    if (x>nxmax) nxmax=x;
                    if (x<nxmin) nxmin=x;
                } 
                if ((x-cx)*(x-cx)+(z-cz)*(z-cz)<=1.) {
                    if (y>nymax) nymax=y;
                    if (y<nymin) nymin=y;
                }
                if ((x-cx)*(x-cx)+(y-cy)*(y-cy)<=1.) {
                    if (z>nzmax) nzmax=z;
                    if (z<nzmin) nzmin=z;
                }
                rsq = (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz);
                if (rsq>rmax*rmax) {
                    rmax = sqrt(rsq);
                    if (rsq>0) {
                        ir = 1./rmax;
                        ax = (x-cx)*ir;
                        ay = (y-cy)*ir;
                        az = (z-cz)*ir;
                    }
                }
                if (rsq<rmin*rmin) {
                    rmin = sqrt(rsq);
                    if (rsq>0) {
                        ir = 1./rmin;
                        bx = (x-cx)*ir;
                        by = (y-cy)*ir;
                        bz = (z-cz)*ir;
                    }
                }
            }

            if (nlist>0) {
                fprintf(ofile,"%f %f %f  %f %f %f  %f %f %f  %f %f %f %f  %f %f %f %f\n",cx,cy,cz,exmax-exmin,eymax-eymin,ezmax-ezmin,nxmax-nxmin,nymax-nymin,nzmax-nzmin,rmax,ax,ay,az,rmin,bx,by,bz);
            } else {
                fprintf(ofile,"0. 0. 0.  0. 0. 0.  0. 0. 0.  0. 0. 0. 0.  0. 0. 0. 0.\n");
            }
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
    return 0;
}
