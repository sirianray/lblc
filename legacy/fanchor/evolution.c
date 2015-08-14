#include "lb.h"
#include "particle.h"

void evol_p(real Q0[Nx][Ny][Nz][3][3])
{
  if(particle_on!=0){
    	p_update();						// update particle position, etc.
		if (p_move!=0) {
		  p_iden_up(Q0);				// if particle has displaced, update link points, etc.
		  if (t_current%2==0)p_link_up();
		}
	}
}

void evol_f(real f1[Nx][Ny][Nz][15], real f2[Nx][Ny][Nz][15], real Q0[Nx][Ny][Nz][3][3])
{
	int i,j,k,ii,jj;
	int ip, jp, kp;
	real Cf1, Cf2, CG1, CG2;
	real temp;

	cal_feqf_new(f_eq, f1, sigma);
	cal_p(p, f1, f_eq);

	streaming(f1,f2);
	if(particle_on!=0)p_bc(f1,f2);
	
	for(i=0;i<Nx;i++){
	    for(j=0;j<Ny;j++){
	    	for(k=0;k<Nz;k++){
				for(ii=0;ii<15;ii++) {
					Cf[i][j][k][ii] = - itau_f * (f1[i][j][k][ii] - f_eq[i][j][k][ii]) + p[i][j][k][ii];
					ip = i + e[ii][0];
					jp = j + e[ii][1];
					kp = k + e[ii][2];					
					if ( ((kp>=0 && kp <Nz) || wall_on==0) && (particle_on==0 || particle_on !=0 && link_point2[i][j][k][ii]==0) ) {
						if (ip >= Nx) {
							ip -= Nx;
						} else if (ip<0) {
							ip += Nx;
						}
						if (jp >= Ny) {
							jp -= Ny;
						} else if (jp<0) {
							jp += Ny;
						}
						jj = ii;
					}
					else {
						ip = i;
						jp = j;
						kp = k;
						jj = bounce(ii);
					}
					f2[ip][jp][kp][jj] += dt * Cf[i][j][k][ii] ;
					
				}
			}
		}
	}

	cal_rho(Rho,f2);
	cal_u(u,f2);
	cal_sigma(Q0);
	cal_feqf_new(f_eq,f2,sigma);
	cal_p(p,f2,f_eq);

	for(k=0;k<Nz;k++){
		for(i=0;i<Nx;i++){
		    for(j=0;j<Ny;j++){
				for(ii=0;ii<15;ii++) {
					ip = i + e[ii][0];
					jp = j + e[ii][1];
					kp = k + e[ii][2];
					if ( ((kp>=0 && kp <Nz) || wall_on==0) && (particle_on==0 || particle_on !=0 && link_point2[i][j][k][ii]==0) ) {
						if (ip >= Nx) {
							ip -= Nx;
						} else if (ip<0) {
							ip += Nx;
						}
						if (jp >= Ny) {
							jp -= Ny;
						} else if (jp<0) {
							jp += Ny;
						}
						jj = ii;
					}
					else {
						ip = i;
						jp = j;
						kp = k;
						jj = bounce(ii);
					}
					Cf1 = Cf[i][j][k][ii];
					Cf2 = - itau_f * (f2[ip][jp][kp][jj] - f_eq[ip][jp][kp][jj]) + p[ip][jp][kp][jj];
					f1[i][j][k][ii] = f1[i][j][k][ii] + dt * 0.5 * (Cf1 + Cf2);
				}
			}
		}
	}


	streaming(f1,f2);
	if(particle_on!=0)p_bc(f1,f2);

	cal_rho(Rho,f2);
	cal_u(u, f2);
//	cal_W(u);
//	cal_sigma();
}


void streaming(real f1[Nx][Ny][Nz][15], real f2[Nx][Ny][Nz][15])
{
	int  i, j, k, ii, iip, jj, kk, in, jn, kn;
	real Bx2;
	real dx, dy, dz, dr2i, rrad;
	int  ipar;
	

	for (ii=0; ii<15; ii++) {
		
		if(ii<7) {
			Bx2 = two3rd;
		} else {
			Bx2 = one12th;
		}
		iip = bounce(ii);
		
		for (k=0; k<Nz; k++) {
			for (i=0; i<Nx; i++) {
				for (j=0; j<Ny; j++) {

						in = i + e[ii][0];
						jn = j + e[ii][1];
						kn = k + e[ii][2];
						if(in==-1)in=Nx-1;
						if(in==Nx)in=0;
						if(jn==-1)jn=Ny-1;
						if(jn==Ny)jn=0;
						
						if (kn>=0 && kn<Nz) {
							if (particle_on==0 || link_point2[i][j][k][ii]==0){
								f2[in][jn][kn][ii]=f1[i][j][k][ii];
							}
						}
						else {
							if (wall_on==0) {
								if(kn==-1)kn=Nz-1;
								if(kn==Nz)kn=0;
								f2[in][jn][kn][ii]=f1[i][j][k][ii];
							} 
							else {
								f2[i][j][k][iip]=f1[i][j][k][ii];
								if (k==0) {
									f2[i][j][k][iip] += -Bx2*Rho[i][j][k]*(uy_bottom*(double)e[ii][1]);
								} else if(k==Nz-1) {
									f2[i][j][k][iip] += -Bx2*Rho[i][j][k]*(uy_top*(double)e[ii][1]);
								}
							} 
						}
				
				}
			}
		}
	}
	
}


/* 
//if the wall's normal is in x-direction
void streaming(real f1[Nx][Ny][Nz][15], real f2[Nx][Ny][Nz][15])
{
	int  i, j, k, ii, iip, jj, kk, in, jn, kn;
	int  pbc_x=0;
	real Bx2;
	real dx, dy, dz, dr2i, rrad;
	int  ipar;
	

	for (ii=0; ii<15; ii++) {
		
		if(ii<7) {
			Bx2 = two3rd;
		} else {
			Bx2 = one12th;
		}
		iip = bounce(ii);
		
		for (k=0; k<Nz; k++) {
			for (i=0; i<Nx; i++) {
				for (j=0; j<Ny; j++) {

						in = i + e[ii][0];
						jn = j + e[ii][1];
						kn = k + e[ii][2];
//						if(in==-1)in=Nx-1;
//						if(in==Nx)in=0;
						if(jn==-1)jn=Ny-1;
						if(jn==Ny)jn=0;
						if(kn==-1)kn=Nz-1;
						if(kn==Nz)kn=0;
						
						if (in>=0 && in<Nx) {
							if (particle_on==0 || link_point2[i][j][k][ii]==0){
								f2[in][jn][kn][ii]=f1[i][j][k][ii];
							}
						}
						else {
							if (pbc_x==1) {
								if(in==-1)in=Nx-1;
								if(in==Nx)in=0;
								f2[in][jn][kn][ii]=f1[i][j][k][ii];
							} 
							else {
								f2[i][j][k][iip]=f1[i][j][k][ii];
								if (i==0) {
									f2[i][j][k][iip] += -Bx2*Rho[i][j][k]*(uy*(double)e[ii][1]);
								}
							} 
						}
				
				}
			}
		}
	}
	
}
*/

// yeomans boundary condtion
/*
void bc(real f1[Nx][Ny][Nz][15])
{
	int i, j;
	real f1w=0, f2w=0, f3=0, f4=0, f6 = 0, f11 = 0, f12 = 0, f13 = 0, f14 = 0;
	real f5=0, f7=0, f8=0, f9=0, f10=0;
	real rho_w;

	for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			  rho_w = Rho[i][j][0];
				f1w = f1[i][j][0][1];
				f2w = f1[i][j][0][2];
				f3  = f1[i][j][0][3];
				f4  = f1[i][j][0][4];
				f6  = f1[i][j][0][6];
				f11 = f1[i][j][0][11];
				f12 = f1[i][j][0][12];
				f13 = f1[i][j][0][13];
				f14 = f1[i][j][0][14];
				
				f1[i][j][0][5 ] = f6;
				f1[i][j][0][7 ] = 0.25*(-f1w-f2w+f3+f4-f11+f12+3*f13+f14+rho_w*uy_bottom);
				f1[i][j][0][8 ] = 0.25*( f1w-f2w-f3+f4+f11-f12+f13+3*f14+rho_w*uy_bottom);
				f1[i][j][0][9 ] = 0.25*( f1w+f2w-f3-f4+3*f11+f12-f13+f14-rho_w*uy_bottom);
				f1[i][j][0][10] = 0.25*(-f1w+f2w+f3-f4+f11+3*f12+f13-f14-rho_w*uy_bottom);
				
				f1w = f1[i][j][Nz-1][1];
				f2w = f1[i][j][Nz-1][2];
				f3  = f1[i][j][Nz-1][3];
				f4  = f1[i][j][Nz-1][4];
				f5  = f1[i][j][Nz-1][5];		//6 ->5
				f7  = f1[i][j][Nz-1][7];		//11->7
				f8  = f1[i][j][Nz-1][8];		//12->8
				f9  = f1[i][j][Nz-1][9];		//13->9
				f10 = f1[i][j][Nz-1][10];	//14->10
				
				
				f1[i][j][Nz-1][6 ] = f5;
				f1[i][j][Nz-1][11] = 0.25*(-f1w-f2w+f3+f4-f7+f8+3*f9+f10+rho_w*uy_top);
				f1[i][j][Nz-1][12] = 0.25*( f1w-f2w-f3+f4+f7-f8+f9+3*f10+rho_w*uy_top);
				f1[i][j][Nz-1][13] = 0.25*( f1w+f2w-f3-f4+3*f7+f8-f9+f10-rho_w*uy_top);
				f1[i][j][Nz-1][14] = 0.25*(-f1w+f2w+f3-f4+f7+3*f8+f9-f10-rho_w*uy_top);			
		}
	}
	
}
*/
