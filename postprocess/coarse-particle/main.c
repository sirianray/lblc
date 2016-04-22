#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "math.h"
#include <unistd.h>

int main(int argc, char *argv[]){

	FILE *grid, *file, *ufile, *rfile, *Qfile, *sfile, *ifile;
	float *g1, *g2, *g3, *v1, *v2, *v3;
	float a[9], w[3], u[3], S, density;
	int *gi;
	int Nx, Ny, Nz, points, lines=0, i, m, info, junk;
	char ch, mfile[256];

	int m_int=1, Nxc, Nyc, Nzc, pointsc, j, k, id;

        grid = fopen("grid.out", "r");
        fscanf(grid,"Nx Ny Nz %d %d %d\n",&Nx,&Ny,&Nz);
        points = Nx*Ny*Nz;

	if (argc>1) m_int = atoi(argv[1]);
	printf("every %d mesh points will be printed\n",m_int);
	Nxc = (Nx-1)/m_int + 1;
	Nyc = (Ny-1)/m_int + 1;
	Nzc = (Nz-1)/m_int + 1;
	pointsc = Nxc*Nyc*Nzc;

	g1 = malloc(points*sizeof(float));
        g2 = malloc(points*sizeof(float));
        g3 = malloc(points*sizeof(float));
        v1 = malloc(points*sizeof(float));
        v2 = malloc(points*sizeof(float));
        v3 = malloc(points*sizeof(float));
	gi = malloc(points*sizeof(int));

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
		fscanf(grid,"%f %f %f %d\n", &g1[i], &g2[i], &g3[i], &gi[i]);
	}
	fclose(grid);

	ufile=fopen("u_3d.out","r");
	rfile=fopen("rho_3d.out","r");
	Qfile=fopen("Q_3d.out","r");
	sfile=fopen("stress_3d.out","r");
	ifile=fopen("type_3d.out","r");
	
	for(m=0;m<lines;m++){
		sprintf(mfile,"movie_%05d.vtk",m);
                file = fopen(mfile,"w");

		fprintf(file,"# vtk DataFile Version 3.0\n");
                fprintf(file,"vtk output\n");
                fprintf(file,"ASCII\n");
                fprintf(file,"DATASET STRUCTURED_GRID\n");

		fprintf(file,"DIMENSIONS\t %d\t %d\t %d\t\n", Nxc, Nyc, Nzc);
                fprintf(file,"POINTS\t %d\t float\n",pointsc);

		for (k=0; k<Nz; k+=m_int) {
			for (j=0; j<Ny; j+=m_int) {
				for (i=0; i<Nx; i+=m_int) {
					id = i + ( j + k*Ny)*Nx;
					fprintf(file,"\t%0.2f\t%0.2f\t%0.2f\n", g1[id], g2[id], g3[id]);
				}
			}
		}

                fprintf(file,"\n");

                fprintf(file,"POINT_DATA\t%d\n",pointsc);

                if(ifile!=NULL){
                        fprintf(file,"SCALARS type int 1\n");
			fprintf(file,"LOOKUP_TABLE default\n");
                        for(id=0;id<points;id++){
				i = id%Nx;
				j =(id/Nx)%Ny;
				k =(id/Nx)/Ny;
                                fscanf(ifile,"%d\n",&gi[id]);
                                if (i%m_int==0 && j%m_int==0 && k%m_int==0) fprintf(file,"\t%d\n",gi[i]);
                        }
			fprintf(file,"\n");
                }

//		fprintf(file,"POINT_DATA\t%d\n",points);
                fprintf(file,"SCALARS S_field float 1\n");
                fprintf(file,"LOOKUP_TABLE default\n");

		for(id=0;id<points;id++){
			i = id%Nx;
			j =(id/Nx)%Ny;
			k =(id/Nx)/Ny;
			fscanf(Qfile,"%f %f %f %f %f\n",&a[0],&a[3],&a[6],&a[4],&a[7]);
			if (i%m_int==0 && j%m_int==0 && k%m_int==0) {
				a[8] = -1*(a[0] + a[4]);
				info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, w);

				if(info>0) {
					printf("failure in eigenvalue routine\n");
					exit(1);
				}

				S = 0.5*3*w[2];
				v1[id] = a[6];
				v2[id] = a[7];
				v3[id] = a[8];

				fprintf(file,"\t%lf\n",S);
			}
		}
		
		fprintf(file,"\n");
                fprintf(file,"VECTORS directors float\n");

		for (k=0; k<Nz; k+=m_int) {
			for (j=0; j<Ny; j+=m_int) {
				for (i=0; i<Nx; i+=m_int) {
					id = i + ( j + k*Ny)*Nx;
					fprintf(file,"\t%f\t%f\t%f\n",v1[id],v2[id],v3[id]);
				}
			}
		}

		fprintf(file,"\n");
                fprintf(file,"VECTORS velocity float\n");

		for(id=0;id<points;id++){
                	fscanf(ufile,"%f %f %f\n",&u[0],&u[1],&u[2]);
			i = id%Nx;
			j =(id/Nx)%Ny;
			k =(id/Nx)/Ny;
                        if (i%m_int==0 && j%m_int==0 && k%m_int==0) fprintf(file,"\t%f\t%f\t%f\n",u[0],u[1],u[2]);
                }

		if(sfile!=NULL){
			fprintf(file,"\n");
                        fprintf(file,"VECTORS stress float\n");
			for(id=0;id<points;id++){
                        	fscanf(sfile,"%f %f %f\n",&u[0],&u[1],&u[2]);
				i = id%Nx;
				j =(id/Nx)%Ny;
				k =(id/Nx)/Ny;
				if (i%m_int==0 && j%m_int==0 && k%m_int==0) fprintf(file,"\t%f\t%f\t%f\n",u[0],u[1],u[2]);
                	}
		}

		if(rfile!=NULL){
			fprintf(file,"\n");
                        fprintf(file,"SCALARS density float\n");
                        fprintf(file,"LOOKUP_TABLE default\n");
			for(id=0;id<points;id++){
                        	fscanf(rfile,"%f\n",&density);
				i = id%Nx;
				j =(id/Nx)%Ny;
				k =(id/Nx)/Ny;
				if (i%m_int==0 && j%m_int==0 && k%m_int==0) fprintf(file,"\t%e\n",density);
                        }
		}

		printf("printed file %d\n",m);
                fclose(file);
	}

	fclose(ufile);
	if(rfile!=NULL)fclose(rfile);
	fclose(Qfile);
	if(sfile!=NULL)fclose(sfile);
	if(ifile!=NULL)fclose(ifile);

	free(g1);
        free(g2);
        free(g3);
        free(v1);
        free(v2);
        free(v3);
	free(gi);

	printf("done!\n");
   
  return 0;
}
