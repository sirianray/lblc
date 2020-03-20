/*
 *  func.c
 *  
 *
 *  Created by Sirius on 2/21/15.
 *  Copyright 2015 home. All rights reserved.
 *
 */


#include "main.h"

double getS(double U)
{
	double S=0;
	S = 0.25+0.75*sqrt(1.0-eight3rd/U);
	return S;
}

double getF0(double U)
{
        double S, p[5]={0.}, trqq, qqq, eld;

        S = getS(U);
        p[0] = S * two3rd;
        p[3] =-S * third;

        trqq = trQQ(p);
        qqq = QQQ(p);
        eld = A_ldg * ( 0.5 * (1.0-third*U)*trqq - third*U*qqq + 0.25*U*trqq*trqq );

        return eld;
}

double OneThirdDelta(int i)
{
	if(i==0 || i==3 || i==5) {
		return third;
	}
	return 0.0;
}

double One3rdDelta(int i, int j)
{
	if(i==j) {
		return third;
	}
	return 0.0;
}

double Delta(int i, int j)
{
	if(i==j){
		return 1.0;
	}
	return 0.0;
}

double trQQ(real *q)
{
	double q0=*q, q1=*(q+1), q2=*(q+2), q3=*(q+3), q4=*(q+4);
	
	return 2.0*(q0*q0+q1*q1+q2*q2+q3*q3+q4*q4+q0*q3);
}

double QQ(real *q, int k)
{
	int i, j, ii;
	double Qm[3][3], sum=0;
	
	Qm[0][0] = *q;
	Qm[0][1] = *(q+1);
	Qm[0][2] = *(q+2);
	Qm[1][1] = *(q+3);
	Qm[1][2] = *(q+4);
	Qm[2][2] =-(*q)-(*(q+3));
	Qm[1][0] = Qm[0][1];
	Qm[2][0] = Qm[0][2];
	Qm[2][1] = Qm[1][2];
	
	switch (k) {
		case 0:
			i = 0;
			j = 0;
			break;
		case 1:
			i = 0;
			j = 1;
			break;
		case 2:
			i = 0;
			j = 2;
			break;
		case 3:
			i = 1;
			j = 1;
			break;
		case 4:
			i = 1;
			j = 2;
			break;
		case 5:
			i = 2;
			j = 2;
			break;
		default:
			printf("error in QQ calculation\n");
			exit(-1);
			break;
	}
	
	for (ii=0; ii<3; ii++) {
		sum += Qm[i][ii] * Qm[ii][j];
	}
	return sum;
}

double QQQ(real *q)
{
	int ii, jj, kk;
	double Qm[3][3], qqq=0;
	
	Qm[0][0] = *q;
	Qm[0][1] = *(q+1);
	Qm[1][0] = *(q+1);
	Qm[0][2] = *(q+2);
	Qm[2][0] = *(q+2);
	Qm[1][1] = *(q+3);
	Qm[1][2] = *(q+4);
	Qm[2][1] = *(q+4);
	Qm[2][2] =-(*q)-(*(q+3));
	
	for (ii=0; ii<3; ii++) {
		for (jj=0; jj<3; jj++) {
			for (kk=0; kk<3; kk++) {
				qqq += Qm[ii][jj]*Qm[jj][kk]*Qm[kk][ii];
			}
		}
	}
	return qqq;
}

void output1(int new, char d, int cx, int cy)
{
	int i, j, k, ii, xlo=0, xhi=Nx-1, ylo=0, yhi=Ny-1, zlo=0, zhi=Nz-1, id;
	FILE *fout[9];

	if (myid == root) {
		if(new==1){
			if (Q_on!=0) {
				fout[0]=fopen("q1.out","w");
				fout[1]=fopen("q2.out","w");
				fout[2]=fopen("q3.out","w");
				fout[3]=fopen("q4.out","w");
				fout[4]=fopen("q5.out","w");
			}
			if (flow_on!=0) {
				fout[5]=fopen("ux.out","w");
				fout[6]=fopen("uy.out","w");
				fout[7]=fopen("uz.out","w");
				fout[8]=fopen("rho.out","w");
			}
		} else {
			if (Q_on!=0) {
				fout[0]=fopen("q1.out","a");
				fout[1]=fopen("q2.out","a");
				fout[2]=fopen("q3.out","a");
				fout[3]=fopen("q4.out","a");
				fout[4]=fopen("q5.out","a");
			}
			if (flow_on!=0) {
				fout[5]=fopen("ux.out","a");
				fout[6]=fopen("uy.out","a");
				fout[7]=fopen("uz.out","a");
				fout[8]=fopen("rho.out","a");
			}
		}
		
		switch (d) {
			case 'x':
				ylo = cx;
				yhi = cx;
				zlo = cy;
				zhi = cy;
				break;
			case 'y':
				xlo = cx;
				xhi = cx;
				zlo = cy;
				zhi = cy;
				break;
			default:
				xlo = cx;
				xhi = cx;
				ylo = cy;
				yhi = cy;
				break;
		}

		for (i=xlo; i<=xhi; i++) {
			for (j=ylo; j<=yhi; j++) {
				for (k=zlo; k<=zhi; k++) {
					id = i + (j+k*Ny)*Nx;
					if (flow_on!=0) {
						for (ii=0; ii<3; ii++) {
							fprintf(fout[5+ii],"%22.15e ",u[id*3+ii]);
						}
						fprintf(fout[8],"%22.15e",Rho[id]);
					}
					if (Q_on!=0) {
						for (ii=0; ii<5; ii++) {
							fprintf(fout[ii],"%22.15e ",Q[id*5+ii]);
						}
					}								
				}
			}
		}

		if (flow_on!=0) {
			for (ii=5; ii<9; ii++) {
				fprintf(fout[ii],"\n");
				fclose(fout[ii]);
			}
		}
		if (Q_on!=0) {
			for (ii=0; ii<5; ii++) {
				fprintf(fout[ii],"\n");
				fclose(fout[ii]);
			}
		}
	}
}

void output3(int new)
{
	FILE *xfile, *yfile, *zfile, *sfile;
	float a[9], w[3], S;
    double nx=1., ny=1., nz=1.;
	int id, info;

	if (myid == root) {
		if (new!=0) {
            xfile= fopen("nx.out","w");
            yfile= fopen("ny.out","w");
            zfile= fopen("nz.out","w");
            sfile= fopen("S.out","w");
		} else {
            xfile= fopen("nx.out","a");
            yfile= fopen("ny.out","a");
            zfile= fopen("nz.out","a");
            sfile= fopen("S.out","a");
		}
		
		for (id=0; id<points*5; id+=5) {
            a[0] = Q[id];
            a[3] = Q[id+1];
            a[6] = Q[id+2];
            a[4] = Q[id+3];
            a[7] = Q[id+4];
            a[8] = -1*(a[0] + a[4]);
            info = LAPACKE_ssyev(LAPACK_COL_MAJOR, 'v', 'u', 3, a, 3, w);

            if(info>0) {
                printf("failure in eigenvalue routine\n");
                exit(1);
            }

            S = 1.5*w[2];
            if (nx*a[6]+ny*a[7]+nz*a[8]>0) {
                nx= a[6];
                ny= a[7];
                nz= a[8];
            } else {
                nx=-a[6];
                ny=-a[7];
                nz=-a[8];
            }
            fprintf(xfile,"%9.7f ",nx);
            fprintf(yfile,"%9.7f ",ny);
            fprintf(zfile,"%9.7f ",nz);
            fprintf(sfile,"%9.7f ",S);
		}
        fprintf(xfile,"\n");
        fprintf(yfile,"\n");
        fprintf(zfile,"\n");
        fprintf(sfile,"\n");
		
		fclose(xfile);
		fclose(yfile);
		fclose(zfile);
		fclose(sfile);
	}
}

void output_time(double t_begin)
{
	double t_end, t_spend;

	t_end = MPI_Wtime();

        if (myid == root) {
                //Calculate time used.
                t_spend = (double)(t_end - t_begin) / 60.0;

                if(t_spend < 60){
                        printf("\nElasped time:    %lf min.\n", t_spend);
                }
                else{  
                        printf("\nElasped time:    %lf h.\n", t_spend / 60.0);
                }
        }
}

void monitor()
{
	double err, k_eng_new=0., rhototal=0, utotal[3]={0.}, k_diff;
	int i, iu;
		
	if(myid==0){
		printf("\nt=%d\n",t_current);
		
		if (Q_on!=0) {
			err = (e_toto-e_tot)/e_tot;
			printf("E_tot=%lf, E_ld=%lf, E_el=%lf, E_ch=%lf, E_sf=%lf\n",e_tot,e_ld,e_el,e_ch,e_sf);
			printf("E_L1=%le, E_L2=%le, E_L3=%le, E_L4=%le\n",e_L1,e_L2,e_L3,e_L4);
			e_toto=e_tot;
			printf("dif_Q^2=%35.30f\n",Q_diff);
			if (Q_diff<Q_tol && Q_diff>-Q_tol ) {
				printf("Q converged\n");
				qconverge=1;
			}
		}
		
		if (flow_on!=0) {
			for (i=0; i<points; i++) {
				iu         = i*3;
				k_eng_new += Rho[i]*(u[iu]*u[iu]+u[iu+1]*u[iu+1]+u[iu+2]*u[iu+2]);
				rhototal  += Rho[i];
				utotal[0] += u[iu];
				utotal[1] += u[iu+1];
				utotal[2] += u[iu+2];
			}
			k_eng_new *= 0.5;
			
			printf("vel=%20.15f %20.15f %20.15f\n",utotal[0], utotal[1], utotal[2]);
			printf("rho=%20.15f\n",rhototal);
			if (k_eng<5e-25 && k_eng>-5e-25) {
				k_diff = k_eng_new-k_eng;
			} else {
				k_diff = (k_eng_new-k_eng)/k_eng;
			}		
			printf("diff_k =%30.25f\n",k_diff);
			if( k_diff<u_tol && k_diff>-u_tol ) {
				printf("u converged\n");
				uconverge=1;
			}		
			k_eng=k_eng_new;
		}
	}
	MPI_Bcast(&uconverge, 1, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Bcast(&qconverge, 1, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
}

void write_restart()
{
	FILE *frestart;

	if (myid == root) {
		frestart=fopen("restart.dat","wb");
		
		fwrite(&Nx,sizeof(int),1,frestart);
		fwrite(&Ny,sizeof(int),1,frestart);
		fwrite(&Nz,sizeof(int),1,frestart);
		fwrite(&wall_x,sizeof(int),1,frestart);
		fwrite(&wall_y,sizeof(int),1,frestart);
		fwrite(&wall_z,sizeof(int),1,frestart);
		fwrite(Q,sizeof(real),5*points,frestart);
		fwrite(Rho,sizeof(real),points,frestart);
		fwrite(u,sizeof(real),3*points,frestart);
		if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) fwrite(Qsurf,sizeof(real),5*nodes,frestart);
	}
}

void read_restart()
{
	int Nx0, Ny0, Nz0, wall_x0, wall_y0, wall_z0;
	FILE *frestart;
	
	if (myid == root) {
		frestart=fopen("restart.dat","rb");
		
		fread(&Nx0,sizeof(int),1,frestart);
		fread(&Ny0,sizeof(int),1,frestart);
		fread(&Nz0,sizeof(int),1,frestart);
		fread(&wall_x0,sizeof(int),1,frestart);
		fread(&wall_y0,sizeof(int),1,frestart);
		fread(&wall_z0,sizeof(int),1,frestart);
		if (Nx0!=Nx || Ny0!=Ny || Nz0!=Nz || wall_x0!=wall_x || wall_y0!=wall_y || wall_z0!=wall_z) {
			printf("restart file with   Nx, Ny, Nz, wall_x, wall_y, wall_z=%d %d %d %d %d %d\n",Nx,Ny,Nz,wall_x, wall_y, wall_z);
			printf("but param file with Nx, Ny, Nz, wall_x, wall_y, wall_z=%d %d %d %d %d %d\n",Nx0,Ny0,Nz0,wall_x0, wall_y0, wall_z0);
			exit(-1);
		}
		
		fread(Q,sizeof(real),points*5,frestart);
		fread(Rho,sizeof(real),points,frestart);
		fread(u,sizeof(real),3*points,frestart);
		if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) fread(Qsurf,sizeof(real),5*nodes,frestart);

	}
	output3(1);
}

void normalize(double *x, double *y, double *z)
{
	double r, nx, ny, nz, ir;
	
	nx = *x;
	ny = *y;
	nz = *z;
	r  = sqrt(nx*nx+ny*ny+nz*nz);
	if (r<1e-22) {
		printf("warning: normal vector too small");
	}
	ir = 1.0/r;
	*x = nx*ir;
	*y = ny*ir;
	*z = nz*ir;
}

int mod(int n, int k)
{
    // if k=5, returns -2 ~ +2
    int m;
    
    m = n%k;
    
    if (m>k/2) m-=k;
    if (m<-k/2) m+=k;
    
    return m;
}
