#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "math.h"
#include <unistd.h>

int main(int argc, char *argv[])
{
    
    FILE *gfile, *qfile, *ofile;
    int Nx, Ny, Nz, mx, my, iarg, xlo=0, xhi, ylo=0, yhi, eof, lines=0;
    int id, idm, idp, info, i, j, im, ip, jm, jp;
    float a[9], w[3];
    double nx0, ny0, nz0, nxm, nym, nzm, nxp, nyp, nzp;
    double nxcx, nxcy, nycx, nycy, nzcx, nzcy, rootx, rooty, rootz;
    double splay, twist, bend, isplay, itwist, ibend2;
    double *nx, *ny, *nz;
    double r2, ir;
    char norm;
        
//  read box size
	gfile  = fopen("grid.out", "r");
	fscanf(gfile,"Nx Ny Nz %d %d %d\n",&Nx,&Ny,&Nz);
    fclose(gfile);
    printf("box size = [%d %d %d]\n",Nx,Ny,Nz);
        
//  read normal direction
    iarg = 1;
    if (argc<2) {
        printf("normal direction needs to be determined\n");
        exit(-1);
    } else {
        norm = *argv[iarg];
        printf("normal direction is %c\n",norm);
        switch (norm) {
            case 'x':
                mx = Ny;
                my = Nz;
                xhi= Ny-1;
                yhi= Nz-1;
                break;
            case 'y':
                mx = Nx;
                my = Nz;
                xhi= Nx-1;
                yhi= Nz-1;
                break;
            default:
                mx = Nx;
                my = Ny;
                xhi= Nx-1;
                yhi= Ny-1;
                break;
        }
        if (argc==6) {
            xlo = atoi(argv[2]);
            xhi = atoi(argv[3]);
            ylo = atoi(argv[4]);
            yhi = atoi(argv[5]);
            if (xlo<0) xlo=0;
            if (ylo<0) ylo=0;
            if (xhi>=mx) xhi=mx-1;
            if (yhi>=my) yhi=my-1;
        }
        printf("box of interest = [%d %d %d %d]\n",xlo,xhi,ylo,yhi);
    }

    printf("mx, my = %d %d\n",mx,my);
    nx    = malloc(mx*my*sizeof(double));
    ny    = malloc(mx*my*sizeof(double));
    nz    = malloc(mx*my*sizeof(double));

    qfile = fopen("Q_2d.out","r");
    ofile = fopen("elastic_2d","w");
    
    eof = 0;
    while (eof>=0) {
		lines++;
		printf("frame %d\n",lines);
        for (id=0; id<mx*my; id++) {
            if (id==mx*my-1) {
                eof = fscanf(qfile,"%f %f %f %f %f \n",&a[0],&a[3],&a[6],&a[4],&a[7]);
            } else eof = fscanf(qfile,"%f %f %f %f %f ",&a[0],&a[3],&a[6],&a[4],&a[7]);
            if (eof>0) {
				a[8] = -1*(a[0] + a[4]);
				info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, w);

				if(info>0) {
					printf("failure in eigenvalue routine\n");
					exit(1);
				}

//				S = 1.5*w[2];
                switch (norm) {
                    case 'x':
                        nx[id] = a[7];
                        ny[id] = a[8];
                        nz[id] = a[6];
                        break;
                    case 'y':
                        nx[id] = a[6];
                        ny[id] = a[8];
                        nz[id] = a[7];
                        break;
                    default:
			            nx[id] = a[6];
			            ny[id] = a[7];
			            nz[id] = a[8];
                        break;
                }

                r2 = nx[id]*nx[id] + ny[id]*ny[id];
                if (r2>0) {
                    ir = sqrt(1./r2);
                    nx[id] *= ir;
                    ny[id] *= ir;
                    nz[id]  = 0.;
                }
            }
        }

        if (eof>0) {
            splay = 0.;
            twist = 0.;
            bend  = 0.;
            for (i=xlo; i<=xhi; i++) {
                for (j=ylo; j<=yhi; j++) {
                    id  = i + j*mx;
                    nx0 = nx[id];
                    ny0 = ny[id];
                    nz0 = nz[id];

                    im  = i - 1;
                    if (im<0) im = mx-1;
                    ip  = i + 1;
                    if (ip>=mx) ip = 0;
                    idm = im+ j*mx;
                    idp = ip+ j*mx;
                    nxm = nx[idm];
                    nym = ny[idm];
                    nzm = nz[idm];
                    if (nx0*nxm+ny0*nym+nz0*nzm<0) {
                        nxm =-nxm;
                        nym =-nym;
                        nzm =-nzm;
                    }
                    nxp = nx[idp];
                    nyp = ny[idp];
                    nzp = nz[idp];
                    if (nx0*nxp+ny0*nyp+nz0*nzp<0) {
                        nxp =-nxp;
                        nyp =-nyp;
                        nzp =-nzp;
                    }
                    nxcx=0.5*(nxp - nxm);
                    nycx=0.5*(nyp - nym);
                    nzcx=0.5*(nzp - nzm);

                    jm = j - 1;
                    if (jm<0) jm = my-1;
                    jp = j + 1;
                    if (jp>=my) jp = 0;
                    idm = i+ jm*mx;
                    idp = i+ jp*mx;
                    nxm = nx[idm];
                    nym = ny[idm];
                    nzm = nz[idm];
                    if (nx0*nxm+ny0*nym+nz0*nzm<0) {
                        nxm =-nxm;
                        nym =-nym;
                        nzm =-nzm;
                    }
                    nxp = nx[idp];
                    nyp = ny[idp];
                    nzp = nz[idp];
                    if (nx0*nxp+ny0*nyp+nz0*nzp<0) {
                        nxp =-nxp;
                        nyp =-nyp;
                        nzp =-nzp;
                    }
                    nxcy=0.5*(nxp - nxm);
                    nycy=0.5*(nyp - nym);
                    nzcy=0.5*(nzp - nzm);

                    isplay = nxcx + nycy;
                    splay += isplay*isplay;

                    rootx  = nzcy;
                    rooty  =-nzcx;
                    rootz  = nycx-nxcy;
                    itwist = nx0*rootx + ny0*rooty + nz0*rootz;
                    twist += itwist*itwist;
                    ibend2 = (ny0*rootz-nz0*rooty)*(ny0*rootz-nz0*rooty);
                    ibend2+= (nz0*rootx-nx0*rootz)*(nz0*rootx-nx0*rootz);
                    ibend2+= (nx0*rooty-ny0*rootx)*(nx0*rooty-ny0*rootx);
                    bend  += ibend2;
                }
            }
            fprintf(ofile,"%f %f %f\n",0.5*splay,0.5*twist,0.5*bend);
        }
    }

    free(nx);
    free(ny);
    free(nz);
    fclose(qfile);
    fclose(ofile);
    return 0;
}
