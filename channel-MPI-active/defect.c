#include "main.h"

void d_allocate()
{
    FILE *dfile;
    int i;
    double dx, dy, dr;

    if (rand_init==-40 && myid==root) {
        dfile = fopen("defect.in","r");
        fscanf(dfile,"ndefect %d xl %lf\n",&ndefect, &d_ixl);
        d_ixl= 1./d_ixl;
        d_x = malloc(ndefect*sizeof(double));
        d_y = malloc(ndefect*sizeof(double));
        d_nx= malloc(ndefect*sizeof(double));
        d_ny= malloc(ndefect*sizeof(double));
        d_ch= malloc(ndefect*sizeof(double));

        for (i=0; i<ndefect; i++) {
            fscanf(dfile,"%lf %lf %lf %lf %lf\n",&d_x[i],&d_y[i],&dx,&dy,&d_ch[i]);
            dr = sqrt(dx*dx+dy*dy);
            if (dr>0 && dr!=1) {
                dx/=dr;
                dy/=dr;
            }
            d_nx[i] = dx;
            d_ny[i] = dy;
//            printf("%d: %f %f %f %f %f\n",i,d_x[i],d_y[i],d_nx[i],d_ny[i],d_ch[i]);
        }

        fclose(dfile);
    }
}


void d_deallocate()
{
    if (rand_init==-40 && myid==root) {
        free(d_x);
        free(d_y);
        free(d_nx);
        free(d_ny);
        free(d_ch);
    }
}
