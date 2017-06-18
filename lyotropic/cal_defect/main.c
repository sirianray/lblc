#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main(){

    FILE *grid, *pfile, *ofile;
    int Nx, Ny, Nz, points, xmin, xmax, ymin, ymax, zmin, zmax, neck;
    int ix, iy, iz, i, flag, frame;
    float imean, x1m, x1p, x2m, x2p, y1m, y1p, y2m, y2p, z1m, z1p, z2m, z2p, xlen, ylen, zlen;
    float *phi;

    grid = fopen("grid.out","r");
    fscanf(grid,"Nx Ny Nz %d %d %d\n",&Nx,&Ny,&Nz);
    points = Nx*Ny*Nz;
    fclose(grid);
    printf("Nx Ny Nz points = %d %d %d %d\n",Nx, Ny, Nz, points);

    phi  = malloc(points*sizeof(float));

    pfile= fopen("phi_3d.out","r");
    ofile= fopen("result.out","w");
  
    frame= 0;
    flag = fscanf(pfile,"%f\n", &phi[0]);
    while (pfile!=NULL && flag>0) {
        xmin = Nx;
        xmax = 0;
        ymin = Ny;
        ymax = 0;
        zmin = Nz;
        zmax = 0;
        for (i=0; i<points; i++) {
            if (phi[i]>=0.5) {
                ix = i%Nx;
                iy = (i/Nx)%Ny;
                iz = (i/Nx)/Ny;
 
                if (xmin > ix) xmin = ix;
                if (xmax < ix) xmax = ix;
                if (ymin > iy) ymin = iy;
                if (ymax < iy) ymax = iy;
                if (zmin > iz) zmin = iz;
                if (zmax < iz) zmax = iz;
            }
            if (i<points-1) flag = fscanf(pfile,"%f\n", &phi[i+1]);
        }
        imean = 0.5*((float)xmin + (float)xmax);
 
        neck = 0;
        x1m  = 0.;
        x1p  = 0.;
        x2m  = 0.;
        x2p  = 0.;
        y1m  = 0.;
        y1p  = 0.;
        y2m  = 0.;
        y2p  = 0.;
        z1m  = 0.;
        z1p  = 0.;
        z2m  = 0.;
        z2p  = 0.;
        for (i=0; i<points; i++) {
            ix = i%Nx;
            iy = (i/Nx)%Ny;
            iz = (i/Nx)/Ny;
            if ((ix==(int)(imean) || ix==(int)(imean+1.)) && phi[i]>=0.5) neck++;
//          determine the accurate x, y, z length
            
            if (ix==xmin   && phi[i]>=0.5) x1m += phi[i];
            if (ix==xmin+1 && phi[i]>=0.5) x1p += phi[i];
            if (ix==xmax   && phi[i]>=0.5) x2m += phi[i];
            if (ix==xmax-1 && phi[i]>=0.5) x2p += phi[i];
            if (iy==ymin   && phi[i]>=0.5) y1m += phi[i];
            if (iy==ymin+1 && phi[i]>=0.5) y1p += phi[i];
            if (iy==ymax   && phi[i]>=0.5) y2m += phi[i];
            if (iy==ymax-1 && phi[i]>=0.5) y2p += phi[i];
            if (iz==zmin   && phi[i]>=0.5) z1m += phi[i];
            if (iz==zmin+1 && phi[i]>=0.5) z1p += phi[i];
            if (iz==zmax   && phi[i]>=0.5) z2m += phi[i];
            if (iz==zmax-1 && phi[i]>=0.5) z2p += phi[i];

        }
        xlen = xmax-xmin + (x2m-0.5)/(x2p-x2m) + (x1m-0.5)/(x1p-x1m);
        ylen = ymax-ymin + (y2m-0.5)/(y2p-y2m) + (y1m-0.5)/(y1p-y1m);
        zlen = zmax-zmin + (z2m-0.5)/(z2p-z2m) + (z1m-0.5)/(z1p-z1m);
 
        printf("frame %d\n",frame);
        fprintf(ofile, "%d %d %d %d  %f %f %f %f\n",xmax-xmin,ymax-ymin,zmax-zmin,neck/2,xlen,ylen,zlen,sqrt((double)neck*2./3.1415926535898));
 
        frame++;
        flag = fscanf(pfile,"%f\n", &phi[0]);
    }

    free(phi);
    fclose(pfile);
    fclose(ofile);
    printf("done!\n");
}
