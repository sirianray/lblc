#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "math.h"
#include <unistd.h>

int main(){

	FILE *grid, *file, *ufile, *rfile, *Qfile, *sfile;
	float *g1, *g2, *g3, *v1, *v2, *v3;
	float a[9], w[3], u[3], S, density;
	int Nx, Ny, Nz, points, lines=0, i, m, info, junk;
	char ch, mfile[256];

        grid = fopen("grid.out", "r");
        fscanf(grid,"Nx Ny Nz %d %d %d\n",&Nx,&Ny,&Nz);
        points = Nx*Ny*Nz;

	g1 = malloc(points*sizeof(float));
        g2 = malloc(points*sizeof(float));
        g3 = malloc(points*sizeof(float));
        v1 = malloc(points*sizeof(float));
        v2 = malloc(points*sizeof(float));
        v3 = malloc(points*sizeof(float));

	printf("READ IN PARAMETERS\n");
	printf("number of points: %d\n",points);

//	Check number of frames
        FILE *check;
        check = fopen("u_3d.out","r");
        while(!feof(check)){
        	ch = fgetc(check);
                if(ch == '\n') lines++;
        }
        fclose(check);

        lines = lines/points;
        printf("Number of movie frames: %d\n", lines);

	for(i=0;i<points;i++){
		fscanf(grid,"%f %f %f %d\n", &g1[i], &g2[i], &g3[i], &junk);
	}
	fclose(grid);

	ufile=fopen("u_3d.out","r");
	rfile=fopen("rho_3d.out","r");
	Qfile=fopen("Q_3d.out","r");
	sfile=fopen("stress_3d.out","r");
	
	for(m=0;m<lines;m++){
		sprintf(mfile,"movie_%05d.vtk",m);
                file = fopen(mfile,"w");

		fprintf(file,"# vtk DataFile Version 3.0\n");
                fprintf(file,"vtk output\n");
                fprintf(file,"ASCII\n");
                fprintf(file,"DATASET STRUCTURED_GRID\n");

		fprintf(file,"DIMENSIONS\t %d\t %d\t %d\t\n",Nx,Ny,Nz);
                fprintf(file,"POINTS\t %d\t float\n",points);

		for(i=0;i<points;i++){
                	fprintf(file,"\t%0.2f\t%0.2f\t%0.2f\n", g1[i], g2[i], g3[i]);
                }

		fprintf(file,"\n");

		fprintf(file,"POINT_DATA\t%d\n",points);
                fprintf(file,"SCALARS S_field float 1\n");
                fprintf(file,"LOOKUP_TABLE default\n");

		for(i=0;i<points;i++){
			fscanf(Qfile,"%f %f %f %f %f\n",&a[0],&a[3],&a[6],&a[4],&a[7]);
                        a[8] = -1*(a[0] + a[4]);
                        info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, w);

                        if(info>0) {
	                        printf("failure in eigenvalue routine\n");
                                exit(1);
                        }

                        S = 0.5*3*w[2];
			v1[i] = a[6];
			v2[i] = a[7];
			v3[i] = a[8];

			fprintf(file,"\t%0.2f\n",S);
		}
		
		fprintf(file,"\n");
                fprintf(file,"VECTORS directors float\n");

		for(i=0;i<points;i++){
			fprintf(file,"\t%0.2f\t%0.2f\t%0.2f\n",v1[i],v2[i],v3[i]);
		}

		fprintf(file,"\n");
                fprintf(file,"VECTORS velocity float\n");

		for(i=0;i<points;i++){
                	fscanf(ufile,"%f %f %f\n",&u[0],&u[1],&u[2]);
                        fprintf(file,"\t%f\t%f\t%f\n",u[0],u[1],u[2]);
                }

		if(sfile!=NULL){
			fprintf(file,"\n");
                        fprintf(file,"VECTORS stress float\n");
			for(i=0;i<points;i++){
                        	fscanf(sfile,"%f %f %f\n",&u[0],&u[1],&u[2]);
                        	fprintf(file,"\t%f\t%f\t%f\n",u[0],u[1],u[2]);
                	}
		}

		if(rfile!=NULL){
			fprintf(file,"\n");
                        fprintf(file,"SCALARS density float\n");
                        fprintf(file,"LOOKUP_TABLE default\n");
			for(i=0;i<points;i++){
                        	fscanf(rfile,"%f\n",&density);
                                fprintf(file,"\t%e\n",density);
                        }
		}

		printf("printed file %d\n",m);
                fclose(file);
	}

	fclose(ufile);
	if(rfile!=NULL)fclose(rfile);
	fclose(Qfile);
	if(sfile!=NULL)fclose(sfile);

	free(g1);
        free(g2);
        free(g3);
        free(v1);
        free(v2);
        free(v3);

	printf("done!\n");
   
  return 0;
}
