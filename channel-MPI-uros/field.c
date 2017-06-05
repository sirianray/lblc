/*
 *  field.c
 *  
 *
 *  Created by Sirius on 12/12/15.
 *  Copyright 2015 home. All rights reserved.
 *
 */

#include "field.h"
#include "main.h"

int f_l, f_d;
double f_ea, f_Em, f_w0, f_lbd, f_betaE, f_betaEsurf;
double f_k, f_zr;
double *f_E, *f_H, *f_T, *f_ht, *f_Tsurf, *f_htsurf;
MPI_Win winf_t, winf_e, winf_tsurf, winf_htsurf;

void f_allocate()
{
    if (field_on==1) {
        if (myid==0) printf("allocating field variables...\n");

        f_H = malloc(5*point*sizeof(double));
        f_ht= malloc(  point*sizeof(double));

//        f_E = malloc(3*point*sizeof(double));
//        f_T = malloc(  point*sizeof(double));
        MPI_Win_allocate_shared(point*sizeof(double),1,MPI_INFO_NULL,shmcomm, &f_T, &winf_t);
        MPI_Win_allocate_shared(3*point*sizeof(double),1,MPI_INFO_NULL,shmcomm, &f_E, &winf_e);
        
	    if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) {
		    MPI_Win_allocate_shared(node*sizeof(real),1,MPI_INFO_NULL, shmcomm, &f_Tsurf, &winf_tsurf);
		    MPI_Win_allocate_shared(node*sizeof(real),1,MPI_INFO_NULL, shmcomm, &f_htsurf, &winf_htsurf);
        }
    }
}


void f_deallocate()
{
    if (field_on==1) {
        if (myid==0) printf("deallocating field variables...\n");

        free(f_H);
        free(f_ht);

//        free(f_E);
//        free(f_T);
        MPI_Win_fence(0, winf_t);
        MPI_Win_free(&winf_t);
        MPI_Win_fence(0, winf_e);
        MPI_Win_free(&winf_e);

	    if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) {
            MPI_Win_fence(0, winf_tsurf);
            MPI_Win_free(&winf_tsurf);
            MPI_Win_fence(0, winf_htsurf);
            MPI_Win_free(&winf_htsurf);
        }
    }
}


void f_init()
{
    FILE *ff;
    int id, i, j, k, jk;
    double dx, dy, dz, th, r2, r, w, iw, iw2, ee, ex, ey, ez, e2;

    if (field_on==1) {
        if (myid==0) printf("initiating field variables...\n");

        ff = fopen("field.in","r");
        fscanf(ff,"eps_a Em %le %le\n", &f_ea, &f_Em);
        fscanf(ff,"betaE betaEsurf %le %le\n", &f_betaE, &f_betaEsurf);
        fscanf(ff,"l p %d %d\n", &f_l, &f_d);
        fscanf(ff,"lambda w0 %le %le", &f_lbd, &f_w0);
        fclose(ff);

        f_k = 6.28318530718 / f_lbd;
        f_zr= 0.5*f_w0*f_w0*f_k;

        for (id=0; id<lpoint; id++) {
            i = (id+myid*point)%Nx;
            jk= (id+myid*point)/Nx;
            j = jk%Ny;
            k = jk/Ny;
            dx= (double)i - 0.5*(double)Nx;
            dy= (double)j - 0.5*(double)Ny;
            dz= (double)k - 0.5*(double)Nz;
            th= atan2(dy,dx);
            r2= dx*dx + dy*dy;
            r = sqrt(r2);
            w = f_w0*sqrt(1+(dz*dz)/(f_zr*f_zr));
            iw= 1./w;
            iw2= iw*iw;

            ee= f_w0*iw*exp(-r2*iw2)*f_Em*cos(f_k*dz);
            ex= ee*cos((double)f_l*th);
            ey= ee*sin((double)f_l*th);
            ez= 0;
            e2= ex*ex+ey*ey+ez*ez;
            
            f_E[id*3]  = ex;
            f_E[id*3+1]= ey;
            f_E[id*3+2]= ez;

            f_H[id*5]  = f_ea*(-ex*ex + third*e2);
            f_H[id*5+1]= f_ea*(-ex*ey);
            f_H[id*5+2]= f_ea*(-ex*ez);
            f_H[id*5+3]= f_ea*(-ey*ey + third*e2);
            f_H[id*5+4]= f_ea*(-ey*ez);

            f_T[id]    = 0.;
        }

	    if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) {
            for (id=0; id<node; id++) {
                f_Tsurf[id] = 0.;
            }
        }
    }
}


void f_eq(double f_dt)
{
    int i=0, flag=1;
    double err;

    if (field_on==1) {
        while (flag==1) {
            f_evolveT(&err, f_dt);
            i++;

            if (i%100==0) {
                err = err/(double)(Nx*Ny*Nz);
                MPI_Bcast(&err, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
                if (myid==0) printf("t=%d Temp err=%e\n",i,err);
                if (err<1e-8 || i>8000) {
                    flag=0;
                }
	            MPI_Barrier(MPI_COMM_WORLD);
            }
        }
        
        if (myid==0) printf("Temp stablised\n");
    }

    f_output3(1);
}


void f_evolveT(double *err, double f_dt)
{
    int id, ii, ip1, ip2, nid[6], ip, ib;
    double ierr, d2xT, d2yT, d2zT, E2;

    ierr= 0.;

    for (id=0; id<lpoint; id++) {
        if (info[id]==-1) {
			for (ii=0; ii<6; ii++) {
				nid[ii] = neighb[id*6+ii];
			}
			
			ip1 = nid[0];
			ip2 = nid[1];
			if (ip1%5==0 && ip2%5==0)      d2xT = f_T[ip2/5] - 2.0 * f_T[id] + f_T[ip1/5];
			else if (ip1%5!=0 && ip2%5==0) d2xT = four3rd*f_T[ip2/5]-4.0*f_T[id]+eight3rd*f_Tsurf[(ip1+1)/5];
			else if (ip1%5==0 && ip2%5!=0) d2xT = eight3rd*f_Tsurf[(ip2+1)/5]-4.0*f_T[id]+four3rd*f_T[ip1/5];
			else if (ip1%5!=0 && ip2%5!=0) d2xT = 4.*(f_Tsurf[(ip2+1)/5] -2.*f_T[id] + f_Tsurf[(ip1+1)/5]);

			ip1 = nid[2];
			ip2 = nid[3];
			if (ip1%5==0 && ip2%5==0)      d2yT = f_T[ip2/5] - 2.0 * f_T[id] + f_T[ip1/5];
			else if (ip1%5!=0 && ip2%5==0) d2yT = four3rd*f_T[ip2/5]-4.0*f_T[id]+eight3rd*f_Tsurf[(ip1+1)/5];
			else if (ip1%5==0 && ip2%5!=0) d2yT = eight3rd*f_Tsurf[(ip2+1)/5]-4.0*f_T[id]+four3rd*f_T[ip1/5];
			else if (ip1%5!=0 && ip2%5!=0) d2yT = 4.*(f_Tsurf[(ip2+1)/5] -2.*f_T[id] + f_Tsurf[(ip1+1)/5]);

			ip1 = nid[4];
			ip2 = nid[5];
			if (ip1%5==0 && ip2%5==0)      d2zT = f_T[ip2/5] - 2.0 * f_T[id] + f_T[ip1/5];
			else if (ip1%5!=0 && ip2%5==0) d2zT = four3rd*f_T[ip2/5]-4.0*f_T[id]+eight3rd*f_Tsurf[(ip1+1)/5];
			else if (ip1%5==0 && ip2%5!=0) d2zT = eight3rd*f_Tsurf[(ip2+1)/5]-4.0*f_T[id]+four3rd*f_T[ip1/5];
			else if (ip1%5!=0 && ip2%5!=0) d2zT = 4.*(f_Tsurf[(ip2+1)/5] -2.*f_T[id] + f_Tsurf[(ip1+1)/5]);

            E2=f_E[id*3]*f_E[id*3]+f_E[id*3+1]*f_E[id*3+1]+f_E[id*3+2]*f_E[id*3+2];
            f_ht[id] = d2xT+d2yT+d2zT + 0.5*f_betaE*E2;
            ierr += f_ht[id]*f_ht[id];
        }
    }

    // caution: only vaid for flat surface
	if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) {
		for (id=0; id<node; id++) {
			if (id+myid*node<nsurf) {
                ib  = id* 6;
                ip  = id*10;
                d2xT= 0.;
                d2yT= 0.;
                d2zT= 0.;
                if (surf[ip+2]!=0) {
                    ip1 = neighbsurf[ib];
                    ip2 = neighbsurf[ib+1];
                    d2xT= four3rd*f_T[ip2/5] + eight3rd*f_Tsurf[id]-4.*f_T[ip1/5];
                }
                if (surf[ip+3]!=0) {
                    ip1 = neighbsurf[ib+2];
                    ip2 = neighbsurf[ib+3];
                    d2yT= four3rd*f_T[ip2/5] + eight3rd*f_Tsurf[id]-4.*f_T[ip1/5];
                }
                if (surf[ip+4]!=0) {
                    ip1 = neighbsurf[ib+4];
                    ip2 = neighbsurf[ib+5];
                    d2zT= four3rd*f_T[ip2/5] + eight3rd*f_Tsurf[id]-4.*f_T[ip1/5];
                    d2zT= eight3rd*f_T[ip1/5] - 4.*f_Tsurf[id] + four3rd*0.;
                }

                ip1 /= 5;
                E2=f_E[ip1*3]*f_E[ip1*3]+f_E[ip1*3+1]*f_E[ip1*3+1]+f_E[ip1*3+2]*f_E[ip1*3+2];
//                E2   = 0.;
//  only evolve temp for bottom wall surface
                if (surf[ip+4]>0) f_htsurf[id] = d2xT + d2yT + d2zT + 0.5*f_betaEsurf*E2;
            }
        }
    }

	MPI_Reduce(&ierr, err, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    for (id=0; id<lpoint; id++) {
        f_T[id] += f_dt * f_ht[id];
    }
	if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) {
		for (id=0; id<node; id++) {
			if (id+myid*node<nsurf) {
                f_Tsurf[id] += f_dt * f_htsurf[id];
            }
        }
    }

	MPI_Win_fence(0,winf_t);
    MPI_Barrier(MPI_COMM_WORLD);
}


void f_output3(int new)
{
    FILE *grid, *tfile, *efile;
    int i, j, k, id;

    if (field_on==1) {
//      for test only        
//        for (id=0; id<lpoint; id++) f_T[id] = f_E[id*3]*f_E[id*3]+f_E[id*3+1]*f_E[id*3+1];
//	    MPI_Win_fence(0,winf_t);

        if (myid == root) {
        	if (new!=0) {
        		grid = fopen("grid.out","w");
        		tfile= fopen("T_3d.out","w");
                efile= fopen("E_3d.out","w");
        		
        		fprintf(grid,"Nx Ny Nz %d %d %d\n",Nx,Ny,Nz);
        		for (k=0; k<Nz; k++) {
        			for (j=0; j<Ny; j++) {
        				for (i=0; i<Nx; i++) {
        					fprintf(grid,"%d %d %d %d\n",i,j,k,1);
        				}
        			}
        		}
        		
        		fclose(grid);
        	} else {
        		tfile= fopen("T_3d.out","a");
                efile= fopen("E_3d.out","a");
        	}
        	
        	for (id=0; id<points; id++) {
                fprintf(tfile,"%e\n",f_T[id]);
                fprintf(efile,"%e %e %e\n",f_E[id*3],f_E[id*3+1],f_E[id*3+2]);
            }
        	
        	fclose(tfile);
        	fclose(efile);
        }
    }
}
