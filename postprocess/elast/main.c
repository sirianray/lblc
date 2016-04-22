/*
 * This code is to calculate elastic energies 
 *
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "math.h"
#include <unistd.h>

int main(int argc, char *argv[]){
    FILE *grid, *file, *Qfile, *efile, *ofile;
    int Nx, Ny, Nz, point, points, eof, info, i;
    int junk, npar=0, ipar, frame=0;
    int *flag;
    float *g1, *g2, *g3, *px, *py, *pz, *pr;
    float a[9], w[3], u[3], S, S0=0.6;
    float isplay, itwist, ibend, isaddle, iela, splay, twist, bend, saddle, ela;
    double dx, dy, dz, dr, fjunk;

    if (argc>=2) S0 = atof(argv[1]);

    file = fopen("param.in","r");
    fscanf(file,"newrun_on %d\n",&junk);
    fscanf(file,"Nx Ny Nz %d %d %d\n",&Nx,&Ny,&Nz);
    fscanf(file,"npar %d\n",&npar);
    fclose(file);

    if (npar>0) {
        file= fopen("ppos.in","r");
        px  = malloc(npar*sizeof(float));
        py  = malloc(npar*sizeof(float));
        pz  = malloc(npar*sizeof(float));
        pr  = malloc(npar*sizeof(float));
        for (ipar=0; ipar<npar; ipar++) {
            fscanf(file,"%f %f %f %f %f %f\n",&px[ipar],&py[ipar],&pz[ipar],&pr[ipar],&fjunk,&fjunk);
        }
        fclose(file);
    }

    grid = fopen("grid.out","r");
    fscanf(grid,"Nx Ny Nz %d %d %d\n",&Nx,&Ny,&Nz);
    points = Nx*Ny*Nz;

	g1  = malloc(points*sizeof(float));
    g2  = malloc(points*sizeof(float));
    g3  = malloc(points*sizeof(float));
    flag= malloc(points*sizeof(int));

	for(i=0; i<points; i++){
		fscanf(grid,"%f %f %f %d\n", &g1[i], &g2[i], &g3[i], &junk);
        flag[i] = -1;
        for (ipar=0; ipar<npar; ipar++) {
            dx = g1[i] - px[ipar];
            dy = g2[i] - py[ipar];
            dz = g3[i] - pz[ipar];
            dr = sqrt(dx*dx+dy*dy+dz*dz);
            if ((pr[ipar]>0 && dr<pr[ipar]) || (pr[ipar]<0 && dr>-pr[ipar])) flag[i] = 0;
        }
	}
	fclose(grid);

	printf("READ IN PARAMETERS\n");
	printf("number of points:     %d\n",points);
    printf("number of confiments: %d\n",npar);
    printf("threshold S0:         %f\n",S0);

    Qfile = fopen("Q_3d.out","r");
    efile = fopen("elastic_3d.out","r");
    ofile = fopen("elast","w");

    if(efile==NULL) {
        printf("elastic_3d not exist, program ends\n");
        fclose(Qfile);
        fclose(efile);
        fclose(ofile);
        exit(0);
    }

    eof   = fscanf(Qfile,"%f %f %f %f %f\n",&a[0],&a[3],&a[6],&a[4],&a[7]);
    eof   = fscanf(efile,"%f %f %f %f\n",&isplay,&itwist,&ibend,&isaddle);
    point = 0;
    splay = 0;
    twist = 0;
    bend  = 0;
    saddle= 0;
    ela   = 0;
    while(eof>0) {
        if (flag[point]==-1) {
            a[8] = -1*(a[0] + a[4]);
            info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, w);

            if(info>0) {
                printf("failure in eigenvalue routine\n");
                exit(1);
            }

            S = 0.5*3*w[2];
            if (S>=S0) {
                splay += isplay;
                twist += itwist;
                bend  += ibend;
                saddle+= isaddle;
                ela   += isplay+itwist+ibend+isaddle;
            }
        }

        point++;
        eof  = fscanf(Qfile,"%f %f %f %f %f\n",&a[0],&a[3],&a[6],&a[4],&a[7]);
        eof  = fscanf(efile,"%f %f %f %f\n",&isplay,&itwist,&ibend,&isaddle);
        if (point==points) {
            point=0;
            printf("frame %d\n",frame);
            fprintf(ofile,"%f %f %f %f %f\n",splay,twist,bend,saddle,ela);
            splay = 0;
            twist = 0;
            bend  = 0;
            saddle= 0;
            ela   = 0;
            frame++;
        }
    }

    fclose(Qfile);
    fclose(efile);
    fclose(ofile);
    free(g1);
    free(g2);
    free(g3);
    free(flag);
    free(px);
    free(py);
    free(pz);
    free(pr);
    return 0;
}
