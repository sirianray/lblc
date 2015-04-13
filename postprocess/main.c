#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "math.h"
#include <unistd.h>

int main(void){

		  int points=0, i=0, junk=0;
		  int Nx, Ny, Nz, info, m=0, lines=0, n_part;
		  char mfile[100], output[100], ch;
		  float w[3],a[9],dir[3], x, y, z, S=0, R=0, x_part=0, y_part=0, z_part=0;

		  //Eigen routine variables
		  for(i=0;i<9;i++){
					 a[i] = 0;
		  }
		  for(i=0;i<3;i++){
					 w[i] = 0;
					 dir[i] = 0;
		  }

		  float Qmat[3][3];

		  //Read System Parameters
		  FILE *grid;
		  grid = fopen("grid", "r");
		  fscanf(grid,"Nx Ny Nz %d %d %d\n",&Nx,&Ny,&Nz);
		  points = Nx*Ny*Nz;
		  rewind(grid);
		  
		  printf("READ IN PARAMETERS\n");

		  //Check number of frames
		  FILE *check;
		  check = fopen("u_3d","r");
		  while(!feof(check)){
					 ch = fgetc(check);
					 if(ch == '\n') lines++;
		  }
		  fclose(check);

		  lines = lines/points;
		  printf("Number of movie frames: %d\n", lines);

		  //Print out vtk headers
		  FILE *file;
		  FILE *partout;
		  for(m=0;m<lines;m++){

					 sprintf(mfile,"movie_%05d.vtk",m);
					 file = fopen(mfile,"w");

					//Print sphere locations
                                         FILE *sphere;
                                         sphere = fopen("trj.part", "r");

					 fprintf(file,"# vtk DataFile Version 3.0\n");
					 fprintf(file,"vtk output\n");
					 fprintf(file,"ASCII\n");
					 fprintf(file,"DATASET STRUCTURED_GRID\n"); 

					 if(sphere!=NULL){
						sprintf(mfile,"part_%05d.vtk",m);
                                         	partout = fopen(mfile,"w");
						fprintf(partout,"# vtk DataFile Version 2.0\n");
					 	fprintf(partout,"Unstructured grid legacy vtk file with point scalar data\n");
					 	fprintf(partout,"ASCII\n");
					}

					 //Grid Points
					 fprintf(file,"DIMENSIONS\t %d\t %d\t %d\t\n",Nx,Ny,Nz);
					 fprintf(file,"POINTS\t %d\t float\n",points);

					 fscanf(grid,"Nx Ny Nz %d %d %d\n", &junk, &junk, &junk); //Ignore first row
					 for(i=0;i<points;i++){
								fscanf(grid,"%f %f %f %d\n", &x, &y, &z, &junk);
								fprintf(file,"\t%0.2f\t%0.2f\t%0.2f\n", x, y, z);
					 }
					 rewind(grid);

					 fprintf(file,"\n");

					 if(sphere!=NULL)fscanf(sphere,"%d\n",&n_part);

					 if(m>0 && sphere!=NULL){
								for(i=0;i<m;i++){
										  fscanf(sphere,"%f %f %f %f\n",&junk,&junk,&junk,&junk);
								}
					 }

					 if(sphere!=NULL){
						fscanf(sphere,"%f %f %f %f\n",&x_part,&y_part,&z_part,&R);

					 	fprintf(partout,"\nDATASET UNSTRUCTURED_GRID\n");
					 	fprintf(partout,"POINTS 2 double\n");

					 	fclose(sphere);

						 fprintf(partout,"%f %f %f\n",0.0,0.0,0.0);
						 fprintf(partout,"%f %f %f\n",x_part,y_part,z_part);
						 fprintf(partout,"\n");

						 fprintf(partout,"\nPOINT_DATA 2\n");
						 fprintf(partout,"SCALARS radii double\n");
						 fprintf(partout,"LOOKUP_TABLE default\n");
						 fprintf(partout,"0\n");
						 fprintf(partout,"%f\n",R);
						 fprintf(partout,"\n");
					}
					 fprintf(file,"POINT_DATA\t%d\n",points);
					 fprintf(file,"SCALARS S_field float 1\n");
					 fprintf(file,"LOOKUP_TABLE default\n");

					 //Eigen vector and value
					 FILE *Qfile;

					 Qfile = fopen("Q_3d","r");

					 //Skip unnecessary lines
					 if(m>0){
								for(i=0;i<points*m;i++){
										  fscanf(Qfile,"%f %f %f %f %f\n",&junk,&junk,&junk,&junk,&junk);
								}
					 }

					 for(i=0;i<points;i++){
S = 1.0;

								//fscanf(grid,"%f %f %f %d\n", &x, &y, &z, &junk);

								fscanf(Qfile,"%f %f %f %f %f\n",&a[0],&a[3],&a[6],&a[4],&a[7]);
								a[8] = -1*(a[0] + a[4]);

								info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, w);

								if(info>0) {
										  printf("failure in eigenvalue routine\n");
										  exit(1);
								}

								S = 0.5*3*w[2];

								fprintf(file,"\t%0.2f\n",S);
					 }
					 fclose(Qfile);

					 //DIRECTOR FIELD
					 for(i=0;i<9;i++){
								a[i] = 0;
					 }
					 for(i=0;i<3;i++){
								w[i] = 0;
								dir[i] = 0;
					 }

					 fprintf(file,"\n");
					 fprintf(file,"VECTORS directors float\n");

					 Qfile = fopen("Q_3d","r");

					 //Skip unnecessary lines
					 if(m>0){
								for(i=0;i<points*m;i++){
										  fscanf(Qfile,"%f %f %f %f %f\n",&junk,&junk,&junk,&junk,&junk);
								}
					 }

					 for(i=0;i<points;i++){
								dir[0] = 0;
								dir[1] = 0;
								dir[2] = 0;

								fscanf(Qfile,"%f %f %f %f %f\n",&a[0],&a[3],&a[6],&a[4],&a[7]);
								a[8] = -1*(a[0] + a[4]);

								info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, w);

								if(info>0) {
										  printf("failure in eigenvalue routine\n");
										  exit(1);
								}

								dir[0] = a[6];
								dir[1] = a[7];
								dir[2] = a[8];

								fprintf(file,"\t%0.2f\t%0.2f\t%0.2f\n",dir[0],dir[1],dir[2]);
					 }
					 fclose(Qfile);

					 //	VELCOCITY FIELD
					 float u[3];
					 int j;
					 for(i=0;i<3;i++) u[i] = 0;

					 FILE *ufile;
					 ufile = fopen("u_3d","r");

					 float u_mag=0, temp=0;

					 if( access("u_3d", R_OK) !=-1 ){

								fprintf(file,"\n");
								fprintf(file,"VECTORS velocity float\n");

								if(m>0){
										  for(i=0;i<points*m;i++){
													 fscanf(ufile,"%f %f %f\n",&junk,&junk,&junk);
										  }
								}


								for(i=0;i<points;i++){

										  for(j=0;j<3;j++) u[j] = 0;

										  fscanf(ufile,"%f %f %f\n",&u[0],&u[1],&u[2]);

										  u_mag = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);

										  if(u_mag > temp) temp = u_mag;
								}

								rewind(ufile);

								if(m>0){
										  for(i=0;i<points*m;i++){
													 fscanf(ufile,"%f %f %f\n",&junk,&junk,&junk);
										  }
								}

								for(i=0;i<points;i++){

										  for(j=0;j<3;j++) u[j] = 0;

										  fscanf(ufile,"%f %f %f\n",&u[0],&u[1],&u[2]);

										  //for(j=0;j<3;j++) u[j] /= temp;

										  fprintf(file,"\t%f\t%f\t%f\n",u[0],u[1],u[2]);

								}

								fclose(ufile);
					 }
					 else if (m==0) printf("Velocity file DNE\n");

					 //STRESS FIELD

					 ufile = fopen("stress_3d","r");
					 if( access("stress_3d", F_OK) !=-1 ){

								fprintf(file,"\n");

								fprintf(file,"VECTORS stress float\n");

								u_mag=0, temp=0;

								if(m>0){
										  for(i=0;i<points*m;i++){
													 fscanf(ufile,"%f %f %f\n",&junk,&junk,&junk);
										  }
								}

								for(i=0;i<points;i++){

										  for(j=0;j<3;j++) u[j] = 0;

										  fscanf(ufile,"%f %f %f\n",&u[0],&u[1],&u[2]);

										  u_mag = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);

										  if(u_mag > temp) temp = u_mag;
								}

								rewind(ufile);

								for(i=0;i<points;i++){

										  for(j=0;j<3;j++) u[j] = 0;

										  fscanf(ufile,"%f %f %f\n",&u[0],&u[1],&u[2]);

										 // for(j=0;j<3;j++) u[j] /= temp;

										  fprintf(file,"\t%f\t%f\t%f\n",u[0],u[1],u[2]);

								}
								fclose(ufile);
					 }
					 else if(m==0) printf("Stress file DNE\n");

					 //Scalar density field
					 ufile = fopen("rho_3d","r");
					 if (access("rho_3d", F_OK) != -1){

					 fprintf(file,"\n");

					 fprintf(file,"SCALARS density float\n");
					 fprintf(file,"LOOKUP_TABLE default\n");

								float density=0;


								if(m>0){
										  for(i=0;i<points*m;i++){
													 fscanf(ufile,"%f %f %f\n",&junk,&junk,&junk);
										  }
								}

								for(i=0;i<points;i++){

										  fscanf(ufile,"%f\n",&density);
										  fprintf(file,"\t%e\n",density);

								}
								fclose(ufile);
					 }
					 else if(m==0) printf("Density file DNE\n");

					 printf("printed file %d\n",m);

					 fclose(file);
					 if(sphere!=NULL)fclose(partout);
		  }


		  printf("DONE\n");

}

