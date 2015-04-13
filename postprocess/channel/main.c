#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "math.h"
#include <unistd.h>

int main(void){

		  FILE *grid, *file, *ufile, *rfile, *Qfile, *sfile;
		  float *g1, *g2, *g3, *v1, *v2, *v3;
		  float a[9], w[3], S, density;
		  double u[3]={0};
		  int Nx, Ny, Nz, points, lines=0, i, m, info, junk, *flag, t, j, k, index;
		  char ch, mfile[256];
		  double ndotv=0; 
		  double mag=0;

		  grid = fopen("grid.out", "r");
		  fscanf(grid,"Nx Ny Nz %d %d %d\n",&Nx,&Ny,&Nz);
		  points = Nx*Ny*Nz;

		  g1 = malloc(points*sizeof(float));
		  g2 = malloc(points*sizeof(float));
		  g3 = malloc(points*sizeof(float));
		  v1 = malloc(points*sizeof(float));
		  v2 = malloc(points*sizeof(float));
		  v3 = malloc(points*sizeof(float));

		  flag = malloc(points*sizeof(int));

		  printf("READ IN PARAMETERS\n");
		  printf("number of points: %d\n",points);

		  //	Check number of frames
		  FILE *check;
		  check = fopen("vel.out","r");
		  if(check!=NULL){
					 while(!feof(check)){
								ch = fgetc(check);
								if(ch == '\n') lines++;
					 }
					 fclose(check);
					 lines = lines/points;
		  } else lines = 1;

		  printf("Number of movie frames: %d\n", lines);

		  for(i=0;i<points;i++){
					 fscanf(grid,"%f %f %f %d\n", &g1[i], &g2[i], &g3[i], &flag[i]);
		  }
		  fclose(grid);

		  ufile=fopen("vel.out","r");
		  rfile=fopen("den.out","r");
		  Qfile=fopen("Q.out","r");
		  sfile=fopen("stress.out","r");

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
								S = 1;

								fscanf(Qfile,"%f %f %f %f %f\n",&a[0],&a[3],&a[6],&a[4],&a[7]);
								a[8] = -1*(a[0] + a[4]);
								info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, w);

								if(info>0) {
										  printf("failure in eigenvalue routine\n");
										  exit(1);
								}

								if(flag[i]==1 || flag[i]==0) S = 0.5*3*w[2];
								//S = 0.5*3*w[2];
								v1[i] = a[6];
								v2[i] = a[7];
								v3[i] = a[8];

								fprintf(file,"\t%0.2f\n",S);
					 }

					 fprintf(file,"\n");
					 fprintf(file,"VECTORS directors float\n");

					 //for(i=0;i<points;i++){
					 for(i=0;i<Nx;i++){
								for(j=0;j<Ny;j++){
										  for(k=0;k<Nz;k++){

													 index = i*Nx+j*Ny+k*Nx*Ny;

													 //if(i==(int)0.5*Nx && j==(int)0.5*Ny) fprintf(file,"\t%0.2f\t%0.2f\t%0.2f\n",v1[index],v2[index],v3[index]);

													 if(flag[i]==1 || flag[i]==0) fprintf(file,"\t%0.2f\t%0.2f\t%0.2f\n",v1[i],v2[i],v3[i]);
													 else fprintf(file,"\t%0.2f\t%0.2f\t%0.2f\n",0.0,0.0,0.0);
										  }
								}
					 }

					 if(ufile!=NULL){
								fprintf(file,"\n");
								fprintf(file,"VECTORS velocity float\n");

								for(i=0;i<points;i++){
										  fscanf(ufile,"%lf %lf %lf\n",&u[0],&u[1],&u[2]);
										  if(flag[i]==1 || flag[i]==0) fprintf(file,"\t%f\t%f\t%f\n",u[0],u[1],u[2]);
										  else fprintf(file,"\t%0.2f\t%0.2f\t%0.2f\n",0.0,0.0,0.0);
								}
					 }

					 if(sfile!=NULL){
								fprintf(file,"\n");
								fprintf(file,"VECTORS stress float\n");
								for(i=0;i<points;i++){
										  fscanf(sfile,"%lf %lf %lf\n",&u[0],&u[1],&u[2]);
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

					 if(ufile!=NULL){
								fprintf(file,"\n");
								fprintf(file,"SCALARS dot float\n");
								fprintf(file,"LOOKUP_TABLE default\n");

								rewind(ufile);

								for(i=0;i<points;i++){
										  ndotv=2;
										  mag = 0;

										  fscanf(ufile,"%lf %lf %lf\n",&u[0],&u[1],&u[2]);

										  if(flag[i]==1 || flag[i]==0){

													 mag = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
													 mag = 1.0/mag;

													 ndotv = mag*(v1[i]*u[0]+v2[i]*u[1]+v3[i]*u[2]);
													 //ndotv= (v1[i]*u[0]+v2[i]*u[1]+v3[i]*u[2]);
													 if(ndotv!=ndotv) printf("ndotv\n");

													 if(ndotv<0) ndotv *= -1;
										  }

										  fprintf(file,"\t%e\n",ndotv);
								}
					 }

					 printf("printed file %d\n",m);
					 fclose(file);
		  }

		  if(ufile!=NULL) fclose(ufile);
		  if(rfile!=NULL) fclose(rfile);
		  if(sfile!=NULL) fclose(sfile);
		  fclose(Qfile);

		  free(g1);
		  free(g2);
		  free(g3);
		  free(v1);
		  free(v2);
		  free(v3);
		  free(flag);

		  printf("Done!\n");
}
