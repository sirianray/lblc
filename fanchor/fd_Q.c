#include "lb.h"
#include "particle.h"

void evol_Q(real Q1[Nx][Ny][Nz][3][3], real Q2[Nx][Ny][Nz][3][3])
{
	int i, j, k, ii, jj, kk;
	real D[3][3], Omega[3][3], S[3][3], trQW, temp;

	cal_dQ(Q1, H, convQ);

	for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			for(k=0;k<Nz;k++){

				trQW=0.0;
				for(ii=0;ii<3;ii++){
					for(jj=0;jj<3;jj++){
						D[ii][jj] = 0.5*(W[i][j][k][ii][jj] + W[i][j][k][jj][ii]);
						Omega[ii][jj] = 0.5*(W[i][j][k][ii][jj] - W[i][j][k][jj][ii]);
						trQW = trQW + Q1[i][j][k][ii][jj]*W[i][j][k][jj][ii];
					 }
		  		}

		  		for(ii=0;ii<3;ii++){
					for(jj=0;jj<3;jj++){
						temp = 0;
						for(kk=0;kk<3;kk++){
							temp += (xi*D[ii][kk] + Omega[ii][kk])*(Q1[i][j][k][kk][jj] + OneThirdDelta(kk,jj));
							temp += (Q1[i][j][k][ii][kk] + OneThirdDelta(ii,kk))*(xi*D[kk][jj] - Omega[kk][jj]);
						}
						temp += - 2.0*xi*(Q1[i][j][k][ii][jj] + OneThirdDelta(ii,jj) )*trQW;
						S[ii][jj] = temp;
					}
		  		}
				
				Q2[i][j][k][0][0] = Q1[i][j][k][0][0] + qdt * (Gamma_rot * H[i][j][k][0][0] + S[0][0] - convQ[i][j][k][0]);
				Q2[i][j][k][0][1] = Q1[i][j][k][0][1] + qdt * (Gamma_rot * H[i][j][k][0][1] + S[0][1] - convQ[i][j][k][1]);
				Q2[i][j][k][0][2] = Q1[i][j][k][0][2] + qdt * (Gamma_rot * H[i][j][k][0][2] + S[0][2] - convQ[i][j][k][2]);
				Q2[i][j][k][1][1] = Q1[i][j][k][1][1] + qdt * (Gamma_rot * H[i][j][k][1][1] + S[1][1] - convQ[i][j][k][3]);
				Q2[i][j][k][1][2] = Q1[i][j][k][1][2] + qdt * (Gamma_rot * H[i][j][k][1][2] + S[1][2] - convQ[i][j][k][4]);
				Q2[i][j][k][1][0] = Q2[i][j][k][0][1];
				Q2[i][j][k][2][0] = Q2[i][j][k][0][2];
				Q2[i][j][k][2][1] = Q2[i][j][k][1][2];
				Q2[i][j][k][2][2] =-Q2[i][j][k][0][0] - Q2[i][j][k][1][1];
				
			}
		}
	}

	//Surface evolution. Assumes inward normal
	int n=0;
	if(W_bot>=0 || W_top>=0 || p_anchor>=0){
		surface_derivative(Q1);
		
		for(i=0;i<Nx;i++){
			for(j=0;j<Ny;j++){
				for(n=0;n<5;n++){
					
					if(W_bot>=0) {
						Qsurf[i][j][0][n]+=-dt*Gamma_rot*(-kappa*dQdxnu[i][j][0][n]+W_bot*(Qsurf[i][j][0][n]-Qsurf0[i][j][0][n]));
					}
					if(W_top>=0) {
						Qsurf[i][j][1][n]+=-dt*Gamma_rot*( kappa*dQdxnu[i][j][1][n]+W_top*(Qsurf[i][j][1][n]-Qsurf0[i][j][1][n]));
					}
				}
			}
		}
		
		if (p_anchor>=0) {
			for (i=0; i<p_nsurf; i++) {
				for(n=0;n<5;n++){
					p_surfQ[i][n] += -dt*Gamma_rot * (-kappa * p_surfdQ[i][n] + p_surfW[i] * (p_surfQ[i][n] - p_surfQ0[i][n]));
				}
			}
		}
	}
	
	
}


void surface_derivative(real Q0[Nx][Ny][Nz][3][3])
{	
	//bottom surface
	int i, j, k, n;
	int i1, j1, k1, i2, j2, k2, ip, jp, kp, nnx, nny, nnz, isurf, jsurf, ii, jj;
	real nx, ny, nz;
	
	real Qvec1[5], Qvec0[5];
	real QvecNx1[5], QvecNx2[5], QvecNy1[5], QvecNy2[5], QvecNz1[5], QvecNz2[5];
	
	for(i=0;i<5;i++){
		Qvec1[i] = 0;
		Qvec0[i] = 0;
		QvecNz1[i] = 0;
		QvecNz2[i] = 0;
	}
	
//	real dx1 = 0.5, dx2 = 1.0, idx = dx1*dx2*(dx1+dx2);
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			
			Qvec0[0] = Q0[i][j][0][0][0];
			Qvec0[1] = Q0[i][j][0][0][1];
			Qvec0[2] = Q0[i][j][0][0][2];
			Qvec0[3] = Q0[i][j][0][1][1];
			Qvec0[4] = Q0[i][j][0][1][2];
			
			Qvec1[0] = Q0[i][j][1][0][0];
			Qvec1[1] = Q0[i][j][1][0][1];
			Qvec1[2] = Q0[i][j][1][0][2];
			Qvec1[3] = Q0[i][j][1][1][1];
			Qvec1[4] = Q0[i][j][1][1][2];
			
			for (n=0; n<5; n++) {
				//Bottom Surface
//dQdxnu[i][j][0][n] = idx*(-dx1*dx1*Qvec1[n] + (dx1+dx2)*(dx1+dx2)*Qvec0[n] - (dx2*dx2+2*dx1*dx2)*Qsurf[i][j][0][n]);
				dQdxnu[i][j][0][n] = -eight3rd*Qsurf[i][j][0][n] + 3.0*Qvec0[n] - third*Qvec1[n];
			}
		}
	}
	
	
//	dx1 = 1.0, dx2 = 0.5, idx = dx1*dx2*(dx1+dx2);
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			
			QvecNz1[0] = Q0[i][j][Nz-1][0][0];
			QvecNz1[1] = Q0[i][j][Nz-1][0][1];
			QvecNz1[2] = Q0[i][j][Nz-1][0][2];
			QvecNz1[3] = Q0[i][j][Nz-1][1][1];
			QvecNz1[4] = Q0[i][j][Nz-1][1][2];
			
			QvecNz2[0] = Q0[i][j][Nz-2][0][0];
			QvecNz2[1] = Q0[i][j][Nz-2][0][1];
			QvecNz2[2] = Q0[i][j][Nz-2][0][2];
			QvecNz2[3] = Q0[i][j][Nz-2][1][1];
			QvecNz2[4] = Q0[i][j][Nz-2][1][2];
			
			for (n=0; n<5; n++) {
				//Top Surface
//dQdxnu[i][j][1][n] = idx*( (dx1*dx1+2*dx1*dx2)*Qsurf[i][j][1][n] - (dx1+dx2)*(dx1+dx2)*QvecNz1[n] + dx2*dx2 * QvecNz2[n]);
				dQdxnu[i][j][1][n] = eight3rd*Qsurf[i][j][1][n] - 3.0*QvecNz1[n] + third*QvecNz2[n];
			}
		}
	}
	
//	finite anchoring particles
	if (p_anchor>=0) {
		for (isurf=0; isurf<p_nsurf; isurf++) {
			ii = p_surf[isurf][3];
			jj = axis_proj(ii);			// dimension
/*			if (ii==1 || ii==3) {
				jj = 0;	// x direction
			} else if (ii==2 || ii==4) {
				jj = 1; // y direction
			} else {
				jj = 2; // z direction
			}*/
			i = (int)p_surf[isurf][0];
			j = (int)p_surf[isurf][1];
			k = (int)p_surf[isurf][2];
			nx= p_surf[isurf][4];
			ny= p_surf[isurf][5];
			nz= p_surf[isurf][6];
			
			if (nx>0) {
				nnx = 1;
			} else if (nx<0) {
				nnx = -1;
			} else {
				nnx = 0;
			}
			if (ny>0) {
				nny = 1;
			} else if (ny<0) {
				nny = -1;
			} else {
				nny = 0;
			}
			if (nz>0) {
				nnz = 1;
			} else if (nz<0) {
				nnz = -1;
			} else {
				nnz = 0;
			}
			ip = i + e[ii][0];
			jp = j + e[ii][1];
			kp = k + e[ii][2];
			i1 = i + nnx;
			j1 = j + nny;
			k1 = k + nnz;
			period_image(&i1,&j1,&k1);
			i2 = i1+ nnx;
			j2 = j1+ nny;
			k2 = k1+ nnz;
			period_image(&i2,&j2,&k2);

			// x direction
			if (jj==0) {
				QvecNx1[0] = Q0[i1][j][k][0][0];
				QvecNx1[1] = Q0[i1][j][k][0][1];
				QvecNx1[2] = Q0[i1][j][k][0][2];
				QvecNx1[3] = Q0[i1][j][k][1][1];
				QvecNx1[4] = Q0[i1][j][k][1][2];
				
				QvecNx2[0] = Q0[i2][j][k][0][0];
				QvecNx2[1] = Q0[i2][j][k][0][1];
				QvecNx2[2] = Q0[i2][j][k][0][2];
				QvecNx2[3] = Q0[i2][j][k][1][1];
				QvecNx2[4] = Q0[i2][j][k][1][2];
				for (n=0; n<5; n++) {
					p_surfdQ[isurf][n] = nnx * nx * (3.0*QvecNx1[n] - third*QvecNx2[n] - eight3rd * p_surfQ[isurf][n]);
				}
			} else if (nnx!=0) {
				jsurf = p_surfi[i1][j][k][jj];
				if (jsurf!=-1) {
					for (n=0; n<5; n++) QvecNx1[n] = p_surfQ[jsurf][n];
				} else {
										
					QvecNx1[0] = 0.5 * (Q0[i1][j][k][0][0] + Q0[i1][jp][kp][0][0]);
					QvecNx1[1] = 0.5 * (Q0[i1][j][k][0][1] + Q0[i1][jp][kp][0][1]);
					QvecNx1[2] = 0.5 * (Q0[i1][j][k][0][2] + Q0[i1][jp][kp][0][2]);
					QvecNx1[3] = 0.5 * (Q0[i1][j][k][1][1] + Q0[i1][jp][kp][1][1]);
					QvecNx1[4] = 0.5 * (Q0[i1][j][k][1][2] + Q0[i1][jp][kp][1][2]);
				}
				
				jsurf = p_surfi[i2][j][k][jj];
				if (jsurf!=-1) {
					for (n=0; n<5; n++) QvecNx2[n] = p_surfQ[jsurf][n];
				} else {
					QvecNx2[0] = 0.5 * (Q0[i2][j][k][0][0] + Q0[i2][jp][kp][0][0]);
					QvecNx2[1] = 0.5 * (Q0[i2][j][k][0][1] + Q0[i2][jp][kp][0][1]);
					QvecNx2[2] = 0.5 * (Q0[i2][j][k][0][2] + Q0[i2][jp][kp][0][2]);
					QvecNx2[3] = 0.5 * (Q0[i2][j][k][1][1] + Q0[i2][jp][kp][1][1]);
					QvecNx2[4] = 0.5 * (Q0[i2][j][k][1][2] + Q0[i2][jp][kp][1][2]);
				}
				for (n=0; n<5; n++) {
					p_surfdQ[isurf][n] = nnx * nx * ( 2.0*QvecNx1[n] - 0.5*QvecNx2[n] - 1.5 * p_surfQ[isurf][n] );
				}
			} else {
				for (n=0; n<5; n++) p_surfdQ[isurf][n] = 0;
			}

			// y direction
			if (jj==1) {
				QvecNy1[0] = Q0[i][j1][k][0][0];
				QvecNy1[1] = Q0[i][j1][k][0][1];
				QvecNy1[2] = Q0[i][j1][k][0][2];
				QvecNy1[3] = Q0[i][j1][k][1][1];
				QvecNy1[4] = Q0[i][j1][k][1][2];
				
				QvecNy2[0] = Q0[i][j2][k][0][0];
				QvecNy2[1] = Q0[i][j2][k][0][1];
				QvecNy2[2] = Q0[i][j2][k][0][2];
				QvecNy2[3] = Q0[i][j2][k][1][1];
				QvecNy2[4] = Q0[i][j2][k][1][2];
				for (n=0; n<5; n++) {
					p_surfdQ[isurf][n] += nny * ny * (3.0*QvecNy1[n] - third*QvecNy2[n] - eight3rd * p_surfQ[isurf][n]);
				}
			} else if (nny!=0) {
				jsurf = p_surfi[i][j1][k][jj];
				if (jsurf!=-1) {
					for (n=0; n<5; n++) QvecNy1[n] = p_surfQ[jsurf][n];
				} else {
					QvecNy1[0] = 0.5 * (Q0[i][j1][k][0][0] + Q0[ip][j1][kp][0][0]);
					QvecNy1[1] = 0.5 * (Q0[i][j1][k][0][1] + Q0[ip][j1][kp][0][1]);
					QvecNy1[2] = 0.5 * (Q0[i][j1][k][0][2] + Q0[ip][j1][kp][0][2]);
					QvecNy1[3] = 0.5 * (Q0[i][j1][k][1][1] + Q0[ip][j1][kp][1][1]);
					QvecNy1[4] = 0.5 * (Q0[i][j1][k][1][2] + Q0[ip][j1][kp][1][2]);
				}
				
				jsurf = p_surfi[i][j2][k][jj];
				if (jsurf!=-1) {
					for (n=0; n<5; n++) QvecNy2[n] = p_surfQ[jsurf][n];
				} else {
					QvecNy2[0] = 0.5 * (Q0[i][j2][k][0][0] + Q0[ip][j2][kp][0][0]);
					QvecNy2[1] = 0.5 * (Q0[i][j2][k][0][1] + Q0[ip][j2][kp][0][1]);
					QvecNy2[2] = 0.5 * (Q0[i][j2][k][0][2] + Q0[ip][j2][kp][0][2]);
					QvecNy2[3] = 0.5 * (Q0[i][j2][k][1][1] + Q0[ip][j2][kp][1][1]);
					QvecNy2[4] = 0.5 * (Q0[i][j2][k][1][2] + Q0[ip][j2][kp][1][2]);
				}
				for (n=0; n<5; n++) {
					p_surfdQ[isurf][n] += nny * ny * ( 2.0*QvecNy1[n] - 0.5*QvecNy2[n] - 1.5 * p_surfQ[isurf][n] );
				}
			}

			
			// z direction
			if (jj==2) {
				QvecNz1[0] = Q0[i][j][k1][0][0];
				QvecNz1[1] = Q0[i][j][k1][0][1];
				QvecNz1[2] = Q0[i][j][k1][0][2];
				QvecNz1[3] = Q0[i][j][k1][1][1];
				QvecNz1[4] = Q0[i][j][k1][1][2];
				
				QvecNz2[0] = Q0[i][j][k2][0][0];
				QvecNz2[1] = Q0[i][j][k2][0][1];
				QvecNz2[2] = Q0[i][j][k2][0][2];
				QvecNz2[3] = Q0[i][j][k2][1][1];
				QvecNz2[4] = Q0[i][j][k2][1][2];
				for (n=0; n<5; n++) {
					p_surfdQ[isurf][n] += nnz * nz * (3.0*QvecNz1[n] - third*QvecNz2[n] - eight3rd * p_surfQ[isurf][n]);
				}
			} else if (nnz!=0) {
				jsurf = p_surfi[i][j][k1][jj];
				if (jsurf!=-1) {
					for (n=0; n<5; n++) QvecNz1[n] = p_surfQ[jsurf][n];
				} else {
					QvecNz1[0] = 0.5 * (Q0[i][j][k1][0][0] + Q0[ip][jp][k1][0][0]);
					QvecNz1[1] = 0.5 * (Q0[i][j][k1][0][1] + Q0[ip][jp][k1][0][1]);
					QvecNz1[2] = 0.5 * (Q0[i][j][k1][0][2] + Q0[ip][jp][k1][0][2]);
					QvecNz1[3] = 0.5 * (Q0[i][j][k1][1][1] + Q0[ip][jp][k1][1][1]);
					QvecNz1[4] = 0.5 * (Q0[i][j][k1][1][2] + Q0[ip][jp][k1][1][2]);
				}
				
				jsurf = p_surfi[i][j][k2][jj];
				if (jsurf!=-1) {
					for (n=0; n<5; n++) QvecNz2[n] = p_surfQ[jsurf][n];
				} else {
					QvecNz2[0] = 0.5 * (Q0[i][j][k2][0][0] + Q0[ip][jp][k2][0][0]);
					QvecNz2[1] = 0.5 * (Q0[i][j][k2][0][1] + Q0[ip][jp][k2][0][1]);
					QvecNz2[2] = 0.5 * (Q0[i][j][k2][0][2] + Q0[ip][jp][k2][0][2]);
					QvecNz2[3] = 0.5 * (Q0[i][j][k2][1][1] + Q0[ip][jp][k2][1][1]);
					QvecNz2[4] = 0.5 * (Q0[i][j][k2][1][2] + Q0[ip][jp][k2][1][2]);
				}
				for (n=0; n<5; n++) {
					p_surfdQ[isurf][n] += nnz * nz * ( 2.0*QvecNz1[n] - 0.5*QvecNz2[n] - 1.5 * p_surfQ[isurf][n] );
				}
			}
			
			
		}
	}
	
}



void cal_dQ(real Q0[Nx][Ny][Nz][3][3], real H0[Nx][Ny][Nz][3][3], real convQ[Nx][Ny][Nz][5])
{

	int i, j, k, ii=0, jj=0, kk=0, isurf, jsurf, scheme=2;
	int ip, im, jp, jm, kp, km;
	int ip2,im2,jp2,jm2,kp2,km2;
	int ipar;
	real t10, t11, t12, t20, t21, t22, t30, t31, t32, t40, t41, t42, t50, t51, t52;
	real temp, Qu[5], Qd[5], Qu2[5], Qd2[5];
	real d2x0, d2x1, d2x2, d2x3, d2x4, d2y0, d2y1, d2y2, d2y3, d2y4, d2z0, d2z1, d2z2, d2z3, d2z4;
	real dx, dy, dz, dr2i;
	real a0, a1, a2, a3, a4;
	real trQ2, trH3rd, d2Q[3][3], dQ2;
	real ux, uy, uz;
	
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for (k=0; k<Nz; k++) {
				a0 = Q0[i][j][k][0][0];
				a1 = Q0[i][j][k][0][1];
				a2 = Q0[i][j][k][0][2];
				a3 = Q0[i][j][k][1][1];
				a4 = Q0[i][j][k][1][2];
                                ux = u[i][j][k][0];
                                uy = u[i][j][k][1];
                                uz = u[i][j][k][2];
				
				im = i-1;
				jm = j-1;
				ip = i+1;
				jp = j+1;
				km = k-1;
				kp = k+1;

				if(i==0) im = Nx-1;
				if(j==0) jm = Ny-1;
				if(k==0) km = Nz-1;

				if(i==Nx-1) ip = 0;
				if(j==Ny-1) jp = 0;
				if(k==Nz-1) kp = 0;

                                im2 = im-1;
                                jm2 = jm-1;
                                ip2 = ip+1;
                                jp2 = jp+1;
                                km2 = km-1;
                                kp2 = kp+1;
                                if(im==0) im2 = Nx-1;
                                if(jm==0) jm2 = Ny-1;
                                if(km==0) km2 = Nz-1;

                                if(ip==Nx-1) ip2 = 0;
                                if(jp==Ny-1) jp2 = 0;
                                if(kp==Nz-1) kp2 = 0;

				if(k==0 && wall_anchoring_on!=0 ){
					Qd[0] = Qsurf[i][j][0][0];
					Qd[1] = Qsurf[i][j][0][1];
					Qd[2] = Qsurf[i][j][0][2];
					Qd[3] = Qsurf[i][j][0][3];
					Qd[4] = Qsurf[i][j][0][4];
				} else {
					Qd[0] = Q0[i][j][km][0][0];
					Qd[1] = Q0[i][j][km][0][1];
					Qd[2] = Q0[i][j][km][0][2];
					Qd[3] = Q0[i][j][km][1][1];
					Qd[4] = Q0[i][j][km][1][2];
				}
				
				if(k==Nz-1 && wall_anchoring_on!=0 ){
					Qu[0] = Qsurf[i][j][1][0];
					Qu[1] = Qsurf[i][j][1][1];
					Qu[2] = Qsurf[i][j][1][2];
					Qu[3] = Qsurf[i][j][1][3];
					Qu[4] = Qsurf[i][j][1][4];
				} else {
					Qu[0] = Q0[i][j][kp][0][0];
					Qu[1] = Q0[i][j][kp][0][1];
					Qu[2] = Q0[i][j][kp][0][2];
					Qu[3] = Q0[i][j][kp][1][1];
					Qu[4] = Q0[i][j][kp][1][2];
				}

                                if(km==0 && wall_anchoring_on!=0 ){
                                        Qd2[0] = Qsurf[i][j][0][0];
                                        Qd2[1] = Qsurf[i][j][0][1];
                                        Qd2[2] = Qsurf[i][j][0][2];
                                        Qd2[3] = Qsurf[i][j][0][3];
                                        Qd2[4] = Qsurf[i][j][0][4];
                                } else {
                                        Qd2[0] = Q0[i][j][km2][0][0];
                                        Qd2[1] = Q0[i][j][km2][0][1];
                                        Qd2[2] = Q0[i][j][km2][0][2];
                                        Qd2[3] = Q0[i][j][km2][1][1];
                                        Qd2[4] = Q0[i][j][km2][1][2];
                                }

                                if(kp==Nz-1 && wall_anchoring_on!=0 ){
                                        Qu2[0] = Qsurf[i][j][1][0];
                                        Qu2[1] = Qsurf[i][j][1][1];
                                        Qu2[2] = Qsurf[i][j][1][2];
                                        Qu2[3] = Qsurf[i][j][1][3];
                                        Qu2[4] = Qsurf[i][j][1][4];
                                } else {
                                        Qu2[0] = Q0[i][j][kp2][0][0];
                                        Qu2[1] = Q0[i][j][kp2][0][1];
                                        Qu2[2] = Q0[i][j][kp2][0][2];
                                        Qu2[3] = Q0[i][j][kp2][1][1];
                                        Qu2[4] = Q0[i][j][kp2][1][2];
                                }
				
				ipar=0;
				if(particle_on && link_point2[i][j][k][1]!=0 && link_point2[i][j][k][3]==0) {
					isurf= p_surfi[i][j][k][0];
					d2x0 = eight3rd*p_surfQ[isurf][0] + four3rd*Q0[im][j][k][0][0] - 4.0 *a0;
					d2x1 = eight3rd*p_surfQ[isurf][1] + four3rd*Q0[im][j][k][0][1] - 4.0 *a1;
					d2x2 = eight3rd*p_surfQ[isurf][2] + four3rd*Q0[im][j][k][0][2] - 4.0 *a2;
					d2x3 = eight3rd*p_surfQ[isurf][3] + four3rd*Q0[im][j][k][1][1] - 4.0 *a3;
					d2x4 = eight3rd*p_surfQ[isurf][4] + four3rd*Q0[im][j][k][1][2] - 4.0 *a4;
					if(scheme==1) {
						t10  = four3rd*p_surfQ[isurf][0] - third*Q0[im][j][k][0][0] - a0;
						t20  = four3rd*p_surfQ[isurf][1] - third*Q0[im][j][k][0][1] - a1;
						t30  = four3rd*p_surfQ[isurf][2] - third*Q0[im][j][k][0][2] - a2;
						t40  = four3rd*p_surfQ[isurf][3] - third*Q0[im][j][k][1][1] - a3;
						t50  = four3rd*p_surfQ[isurf][4] - third*Q0[im][j][k][1][2] - a4;
					} else if(scheme==2) {
						if(ux>0){
							t10 = 1.5 * a0 - 2.0 * Q0[im][j][k][0][0] + 0.5 * Q0[im2][j][k][0][0];
                                                	t20 = 1.5 * a1 - 2.0 * Q0[im][j][k][0][1] + 0.5 * Q0[im2][j][k][0][1];
	                                                t30 = 1.5 * a2 - 2.0 * Q0[im][j][k][0][2] + 0.5 * Q0[im2][j][k][0][2];
        	                                        t40 = 1.5 * a3 - 2.0 * Q0[im][j][k][1][1] + 0.5 * Q0[im2][j][k][1][1];
                	                                t50 = 1.5 * a4 - 2.0 * Q0[im][j][k][1][2] + 0.5 * Q0[im2][j][k][1][2];
						} else if(ux<0){
							t10 = 2.0 * ( p_surfQ[isurf][0] - a0 );
                                        	        t20 = 2.0 * ( p_surfQ[isurf][1] - a1 );
                                                	t30 = 2.0 * ( p_surfQ[isurf][2] - a2 );
	                                                t40 = 2.0 * ( p_surfQ[isurf][3] - a3 );
        	                                        t50 = 2.0 * ( p_surfQ[isurf][4] - a4 );
						}
					}
				} else if(particle_on && link_point2[i][j][k][3]!=0 && link_point2[i][j][k][1]==0) {
					isurf= p_surfi[i-1][j][k][0];
					d2x0 = eight3rd*p_surfQ[isurf][0] + four3rd*Q0[ip][j][k][0][0] - 4.0 *a0;
					d2x1 = eight3rd*p_surfQ[isurf][1] + four3rd*Q0[ip][j][k][0][1] - 4.0 *a1;
					d2x2 = eight3rd*p_surfQ[isurf][2] + four3rd*Q0[ip][j][k][0][2] - 4.0 *a2;
					d2x3 = eight3rd*p_surfQ[isurf][3] + four3rd*Q0[ip][j][k][1][1] - 4.0 *a3;
					d2x4 = eight3rd*p_surfQ[isurf][4] + four3rd*Q0[ip][j][k][1][2] - 4.0 *a4;
					if(scheme==1) {
						t10  =-four3rd*p_surfQ[isurf][0] + third*Q0[ip][j][k][0][0] + a0;
						t20  =-four3rd*p_surfQ[isurf][1] + third*Q0[ip][j][k][0][1] + a1;
						t30  =-four3rd*p_surfQ[isurf][2] + third*Q0[ip][j][k][0][2] + a2;
						t40  =-four3rd*p_surfQ[isurf][3] + third*Q0[ip][j][k][1][1] + a3;
						t50  =-four3rd*p_surfQ[isurf][4] + third*Q0[ip][j][k][1][2] + a4;
					} else if(scheme==2) {
	                                        if(ux<0){
        	                                        t10 =-1.5 * a0 + 2.0 * Q0[ip][j][k][0][0] - 0.5 * Q0[ip2][j][k][0][0];
                	                                t20 =-1.5 * a1 + 2.0 * Q0[ip][j][k][0][1] - 0.5 * Q0[ip2][j][k][0][1];
	                                                t30 =-1.5 * a2 + 2.0 * Q0[ip][j][k][0][2] - 0.5 * Q0[ip2][j][k][0][2];
        	                                        t40 =-1.5 * a3 + 2.0 * Q0[ip][j][k][1][1] - 0.5 * Q0[ip2][j][k][1][1];
                	                                t50 =-1.5 * a4 + 2.0 * Q0[ip][j][k][1][2] - 0.5 * Q0[ip2][j][k][1][2];
                        	                } else if(ux>0){
                                	                t10 = -2.0 * ( p_surfQ[isurf][0] - a0 );
                                        	        t20 = -2.0 * ( p_surfQ[isurf][1] - a1 );
                                                	t30 = -2.0 * ( p_surfQ[isurf][2] - a2 );
	                                                t40 = -2.0 * ( p_surfQ[isurf][3] - a3 );
        	                                        t50 = -2.0 * ( p_surfQ[isurf][4] - a4 );
						}
                                        }
					/*
				} else if (particle_on && link_point2[i][j][k][1]!=0 && link_point2[i][j][k][3]!=0) {
					isurf= p_surfi[i-1][j][k][0];
					jsurf= p_surfi[i  ][j][k][0];
					d2x0 = 4.0 * (p_surfQ[isurf][0] + p_surfQ[jsurf][0] - 2.0 * a0);
					d2x1 = 4.0 * (p_surfQ[isurf][1] + p_surfQ[jsurf][1] - 2.0 * a1);
					d2x2 = 4.0 * (p_surfQ[isurf][2] + p_surfQ[jsurf][2] - 2.0 * a2);
					d2x3 = 4.0 * (p_surfQ[isurf][3] + p_surfQ[jsurf][3] - 2.0 * a3);
					d2x4 = 4.0 * (p_surfQ[isurf][4] + p_surfQ[jsurf][4] - 2.0 * a4);
					t10  = p_surfQ[jsurf][0] - p_surfQ[isurf][0];
					t20  = p_surfQ[jsurf][1] - p_surfQ[isurf][1];
					t30  = p_surfQ[jsurf][2] - p_surfQ[isurf][2];
					t40  = p_surfQ[jsurf][3] - p_surfQ[isurf][3];
					t50  = p_surfQ[jsurf][4] - p_surfQ[isurf][4];
					 */
				} else {
					d2x0 = Q0[ip][j][k][0][0] + Q0[im][j][k][0][0] - 2.0*a0;
					d2x1 = Q0[ip][j][k][0][1] + Q0[im][j][k][0][1] - 2.0*a1;
					d2x2 = Q0[ip][j][k][0][2] + Q0[im][j][k][0][2] - 2.0*a2;
					d2x3 = Q0[ip][j][k][1][1] + Q0[im][j][k][1][1] - 2.0*a3;
					d2x4 = Q0[ip][j][k][1][2] + Q0[im][j][k][1][2] - 2.0*a4;
					if(scheme==1) {
						t10  = 0.5 * (Q0[ip][j][k][0][0]-Q0[im][j][k][0][0]);
						t20  = 0.5 * (Q0[ip][j][k][0][1]-Q0[im][j][k][0][1]);
						t30  = 0.5 * (Q0[ip][j][k][0][2]-Q0[im][j][k][0][2]);
						t40  = 0.5 * (Q0[ip][j][k][1][1]-Q0[im][j][k][1][1]);
						t50  = 0.5 * (Q0[ip][j][k][1][2]-Q0[im][j][k][1][2]);
					} else if(scheme==2) {
	                                        if(ux>0){
        	                                        t10 = 1.5 * a0 - 2.0 * Q0[im][j][k][0][0] + 0.5 * Q0[im2][j][k][0][0];
                	                                t20 = 1.5 * a1 - 2.0 * Q0[im][j][k][0][1] + 0.5 * Q0[im2][j][k][0][1];
                        	                        t30 = 1.5 * a2 - 2.0 * Q0[im][j][k][0][2] + 0.5 * Q0[im2][j][k][0][2];
                                	                t40 = 1.5 * a3 - 2.0 * Q0[im][j][k][1][1] + 0.5 * Q0[im2][j][k][1][1];
                                        	        t50 = 1.5 * a4 - 2.0 * Q0[im][j][k][1][2] + 0.5 * Q0[im2][j][k][1][2];
	                                        } else if(ux<0){
        	                                        t10 =-1.5 * a0 + 2.0 * Q0[ip][j][k][0][0] - 0.5 * Q0[ip2][j][k][0][0];
                	                                t20 =-1.5 * a1 + 2.0 * Q0[ip][j][k][0][1] - 0.5 * Q0[ip2][j][k][0][1];
                        	                        t30 =-1.5 * a2 + 2.0 * Q0[ip][j][k][0][2] - 0.5 * Q0[ip2][j][k][0][2];
                                	                t40 =-1.5 * a3 + 2.0 * Q0[ip][j][k][1][1] - 0.5 * Q0[ip2][j][k][1][1];
                                        	        t50 =-1.5 * a4 + 2.0 * Q0[ip][j][k][1][2] - 0.5 * Q0[ip2][j][k][1][2];
						}
                                        }	
				}
				
				if(particle_on && link_point2[i][j][k][2]!=0 && link_point2[i][j][k][4]==0) {
					isurf= p_surfi[i][j][k][1];
					d2y0 = eight3rd*p_surfQ[isurf][0] + four3rd*Q0[i][jm][k][0][0] - 4.0 *a0;
					d2y1 = eight3rd*p_surfQ[isurf][1] + four3rd*Q0[i][jm][k][0][1] - 4.0 *a1;
					d2y2 = eight3rd*p_surfQ[isurf][2] + four3rd*Q0[i][jm][k][0][2] - 4.0 *a2;
					d2y3 = eight3rd*p_surfQ[isurf][3] + four3rd*Q0[i][jm][k][1][1] - 4.0 *a3;
					d2y4 = eight3rd*p_surfQ[isurf][4] + four3rd*Q0[i][jm][k][1][2] - 4.0 *a4;
					if(scheme==1) {
						t11  = four3rd*p_surfQ[isurf][0] - third*Q0[i][jm][k][0][0] - a0;
						t21  = four3rd*p_surfQ[isurf][1] - third*Q0[i][jm][k][0][1] - a1;
						t31  = four3rd*p_surfQ[isurf][2] - third*Q0[i][jm][k][0][2] - a2;
						t41  = four3rd*p_surfQ[isurf][3] - third*Q0[i][jm][k][1][1] - a3;
						t51  = four3rd*p_surfQ[isurf][4] - third*Q0[i][jm][k][1][2] - a4;
					} else if(scheme==2) {
	                                        if(uy>0){
        	                                        t11 = 1.5 * a0 - 2.0 * Q0[i][jm][k][0][0] + 0.5 * Q0[i][jm2][k][0][0];
                	                                t21 = 1.5 * a1 - 2.0 * Q0[i][jm][k][0][1] + 0.5 * Q0[i][jm2][k][0][1];
                        	                        t31 = 1.5 * a2 - 2.0 * Q0[i][jm][k][0][2] + 0.5 * Q0[i][jm2][k][0][2];
                                	                t41 = 1.5 * a3 - 2.0 * Q0[i][jm][k][1][1] + 0.5 * Q0[i][jm2][k][1][1];
                                        	        t51 = 1.5 * a4 - 2.0 * Q0[i][jm][k][1][2] + 0.5 * Q0[i][jm2][k][1][2];
	                                        } else if(uy<0){
        	                                        t11 = 2.0 * ( p_surfQ[isurf][0] - a0 );
                	                                t21 = 2.0 * ( p_surfQ[isurf][1] - a1 );
                        	                        t31 = 2.0 * ( p_surfQ[isurf][2] - a2 );
                                	                t41 = 2.0 * ( p_surfQ[isurf][3] - a3 );
                                        	        t51 = 2.0 * ( p_surfQ[isurf][4] - a4 );
						}
                                        }
				} else if(particle_on && link_point2[i][j][k][4]!=0 && link_point2[i][j][k][2]==0) {
					isurf= p_surfi[i][j-1][k][1];
					d2y0 = eight3rd*p_surfQ[isurf][0] + four3rd*Q0[i][jp][k][0][0] - 4.0 *a0;
					d2y1 = eight3rd*p_surfQ[isurf][1] + four3rd*Q0[i][jp][k][0][1] - 4.0 *a1;
					d2y2 = eight3rd*p_surfQ[isurf][2] + four3rd*Q0[i][jp][k][0][2] - 4.0 *a2;
					d2y3 = eight3rd*p_surfQ[isurf][3] + four3rd*Q0[i][jp][k][1][1] - 4.0 *a3;
					d2y4 = eight3rd*p_surfQ[isurf][4] + four3rd*Q0[i][jp][k][1][2] - 4.0 *a4;
					if(scheme==1) {
						t11 = -four3rd*p_surfQ[isurf][0] + third*Q0[i][jp][k][0][0] + a0;
						t21 = -four3rd*p_surfQ[isurf][1] + third*Q0[i][jp][k][0][1] + a1;
						t31 = -four3rd*p_surfQ[isurf][2] + third*Q0[i][jp][k][0][2] + a2;
						t41 = -four3rd*p_surfQ[isurf][3] + third*Q0[i][jp][k][1][1] + a3;
						t51 = -four3rd*p_surfQ[isurf][4] + third*Q0[i][jp][k][1][2] + a4;
					} else if(scheme==2) {
	                                        if(uy<0){
        	                                        t11 =-1.5 * a0 + 2.0 * Q0[i][jp][k][0][0] - 0.5 * Q0[i][jp2][k][0][0];
                	                                t21 =-1.5 * a1 + 2.0 * Q0[i][jp][k][0][1] - 0.5 * Q0[i][jp2][k][0][1];
                        	                        t31 =-1.5 * a2 + 2.0 * Q0[i][jp][k][0][2] - 0.5 * Q0[i][jp2][k][0][2];
                                	                t41 =-1.5 * a3 + 2.0 * Q0[i][jp][k][1][1] - 0.5 * Q0[i][jp2][k][1][1];
                                        	        t51 =-1.5 * a4 + 2.0 * Q0[i][jp][k][1][2] - 0.5 * Q0[i][jp2][k][1][2];
                                        	} else if(uy>0){
                                                	t11 =-2.0 * ( p_surfQ[isurf][0] - a0 );
	                                                t21 =-2.0 * ( p_surfQ[isurf][1] - a1 );
        	                                        t31 =-2.0 * ( p_surfQ[isurf][2] - a2 );
                	                                t41 =-2.0 * ( p_surfQ[isurf][3] - a3 );
                        	                        t51 =-2.0 * ( p_surfQ[isurf][4] - a4 );
						}
                                        }
					/*
				} else if (particle_on && link_point2[i][j][k][2]!=0 && link_point2[i][j][k][4]!=0) {
					isurf= p_surfi[i][j-1][k][0];
					jsurf= p_surfi[i][j  ][k][0];
					d2y0 = 4.0 * (p_surfQ[isurf][0] + p_surfQ[jsurf][0] - 2.0 * a0);
					d2y1 = 4.0 * (p_surfQ[isurf][1] + p_surfQ[jsurf][1] - 2.0 * a1);
					d2y2 = 4.0 * (p_surfQ[isurf][2] + p_surfQ[jsurf][2] - 2.0 * a2);
					d2y3 = 4.0 * (p_surfQ[isurf][3] + p_surfQ[jsurf][3] - 2.0 * a3);
					d2y4 = 4.0 * (p_surfQ[isurf][4] + p_surfQ[jsurf][4] - 2.0 * a4);
					t11  = p_surfQ[jsurf][0] - p_surfQ[isurf][0];
					t21  = p_surfQ[jsurf][1] - p_surfQ[isurf][1];
					t31  = p_surfQ[jsurf][2] - p_surfQ[isurf][2];
					t41  = p_surfQ[jsurf][3] - p_surfQ[isurf][3];
					t51  = p_surfQ[jsurf][4] - p_surfQ[isurf][4];	
					 */
				} else {
					d2y0 = Q0[i][jm][k][0][0]+Q0[i][jp][k][0][0]-2.0*a0;
					d2y1 = Q0[i][jm][k][0][1]+Q0[i][jp][k][0][1]-2.0*a1;
					d2y2 = Q0[i][jm][k][0][2]+Q0[i][jp][k][0][2]-2.0*a2;
					d2y3 = Q0[i][jm][k][1][1]+Q0[i][jp][k][1][1]-2.0*a3;
					d2y4 = Q0[i][jm][k][1][2]+Q0[i][jp][k][1][2]-2.0*a4;
					if(scheme==1) {
						t11 = 0.5 * (Q0[i][jp][k][0][0]-Q0[i][jm][k][0][0]);
						t21 = 0.5 * (Q0[i][jp][k][0][1]-Q0[i][jm][k][0][1]);
						t31 = 0.5 * (Q0[i][jp][k][0][2]-Q0[i][jm][k][0][2]);
						t41 = 0.5 * (Q0[i][jp][k][1][1]-Q0[i][jm][k][1][1]);
						t51 = 0.5 * (Q0[i][jp][k][1][2]-Q0[i][jm][k][1][2]);
					} else if(scheme==2) {
	                                        if(uy>0){
        	                                        t11 = 1.5 * a0 - 2.0 * Q0[i][jm][k][0][0] + 0.5 * Q0[i][jm2][k][0][0];
                	                                t21 = 1.5 * a1 - 2.0 * Q0[i][jm][k][0][1] + 0.5 * Q0[i][jm2][k][0][1];
                        	                        t31 = 1.5 * a2 - 2.0 * Q0[i][jm][k][0][2] + 0.5 * Q0[i][jm2][k][0][2];
                                	                t41 = 1.5 * a3 - 2.0 * Q0[i][jm][k][1][1] + 0.5 * Q0[i][jm2][k][1][1];
                                        	        t51 = 1.5 * a4 - 2.0 * Q0[i][jm][k][1][2] + 0.5 * Q0[i][jm2][k][1][2];
                                        	} else if(uy<0){
                                                	t11 =-1.5 * a0 + 2.0 * Q0[i][jp][k][0][0] - 0.5 * Q0[i][jp2][k][0][0];
	                                                t21 =-1.5 * a1 + 2.0 * Q0[i][jp][k][0][1] - 0.5 * Q0[i][jp2][k][0][1];
        	                                        t31 =-1.5 * a2 + 2.0 * Q0[i][jp][k][0][2] - 0.5 * Q0[i][jp2][k][0][2];
                	                                t41 =-1.5 * a3 + 2.0 * Q0[i][jp][k][1][1] - 0.5 * Q0[i][jp2][k][1][1];
                        	                        t51 =-1.5 * a4 + 2.0 * Q0[i][jp][k][1][2] - 0.5 * Q0[i][jp2][k][1][2];
						}
                                        }
				}
				
				if(particle_on && link_point2[i][j][k][5]!=0 && link_point2[i][j][k][6]==0) {
					isurf= p_surfi[i][j][k][2];
					d2z0 = eight3rd*p_surfQ[isurf][0] + four3rd * Qd[0] - 4.0 * a0;
					d2z1 = eight3rd*p_surfQ[isurf][1] + four3rd * Qd[1] - 4.0 * a1;
					d2z2 = eight3rd*p_surfQ[isurf][2] + four3rd * Qd[2] - 4.0 * a2;
					d2z3 = eight3rd*p_surfQ[isurf][3] + four3rd * Qd[3] - 4.0 * a3;
					d2z4 = eight3rd*p_surfQ[isurf][4] + four3rd * Qd[4] - 4.0 * a4;
					if(scheme==1) {
						t12 = four3rd * p_surfQ[isurf][0] - third * Qd[0] - a0;
						t22 = four3rd * p_surfQ[isurf][1] - third * Qd[1] - a1;
						t32 = four3rd * p_surfQ[isurf][2] - third * Qd[2] - a2;
						t42 = four3rd * p_surfQ[isurf][3] - third * Qd[3] - a3;
						t52 = four3rd * p_surfQ[isurf][4] - third * Qd[4] - a4;
					} else if(scheme==2) {
	                                        if(uz>0){
        	                                        t12 = 1.5 * a0 - 2.0 * Q0[i][j][km][0][0] + 0.5 * Q0[i][j][km2][0][0];
                	                                t22 = 1.5 * a1 - 2.0 * Q0[i][j][km][0][1] + 0.5 * Q0[i][j][km2][0][1];
                        	                        t32 = 1.5 * a2 - 2.0 * Q0[i][j][km][0][2] + 0.5 * Q0[i][j][km2][0][2];
                                	                t42 = 1.5 * a3 - 2.0 * Q0[i][j][km][1][1] + 0.5 * Q0[i][j][km2][1][1];
                                        	        t52 = 1.5 * a4 - 2.0 * Q0[i][j][km][1][2] + 0.5 * Q0[i][j][km2][1][2];
	                                        } else if(uz<0){
        	                                        t12 = 2.0 * ( p_surfQ[isurf][0] - a0 );
                	                                t22 = 2.0 * ( p_surfQ[isurf][1] - a1 );
                        	                        t32 = 2.0 * ( p_surfQ[isurf][2] - a2 );
                                	                t42 = 2.0 * ( p_surfQ[isurf][3] - a3 );
                                        	        t52 = 2.0 * ( p_surfQ[isurf][4] - a4 );
						}
                                        }
				} else if(particle_on && link_point2[i][j][k][6]!=0 && link_point2[i][j][k][5]==0) {
					isurf= p_surfi[i][j][k-1][2];
					d2z0 = eight3rd*p_surfQ[isurf][0] + four3rd * Qu[0] - 4.0 * a0;
					d2z1 = eight3rd*p_surfQ[isurf][1] + four3rd * Qu[1] - 4.0 * a1;
					d2z2 = eight3rd*p_surfQ[isurf][2] + four3rd * Qu[2] - 4.0 * a2;
					d2z3 = eight3rd*p_surfQ[isurf][3] + four3rd * Qu[3] - 4.0 * a3;
					d2z4 = eight3rd*p_surfQ[isurf][4] + four3rd * Qu[4] - 4.0 * a4;
					if(scheme==1) {
						t12 =-four3rd * p_surfQ[isurf][0] + third * Qu[0] + a0;
						t22 =-four3rd * p_surfQ[isurf][1] + third * Qu[1] + a1;
						t32 =-four3rd * p_surfQ[isurf][2] + third * Qu[2] + a2;
						t42 =-four3rd * p_surfQ[isurf][3] + third * Qu[3] + a3;
						t52 =-four3rd * p_surfQ[isurf][4] + third * Qu[4] + a4;
					} else if(scheme==2) {
	                                        if(uz<0){
        	                                        t12 =-1.5 * a0 + 2.0 * Q0[i][j][kp][0][0] - 0.5 * Q0[i][j][kp2][0][0];
                	                                t22 =-1.5 * a1 + 2.0 * Q0[i][j][kp][0][1] - 0.5 * Q0[i][j][kp2][0][1];
                        	                        t32 =-1.5 * a2 + 2.0 * Q0[i][j][kp][0][2] - 0.5 * Q0[i][j][kp2][0][2];
                                	                t42 =-1.5 * a3 + 2.0 * Q0[i][j][kp][1][1] - 0.5 * Q0[i][j][kp2][1][1];
                                        	        t52 =-1.5 * a4 + 2.0 * Q0[i][j][kp][1][2] - 0.5 * Q0[i][j][kp2][1][2];
	                                        } else if(uz>0){
        	                                        t12 =-2.0 * ( p_surfQ[isurf][0] - a0 );
                	                                t22 =-2.0 * ( p_surfQ[isurf][1] - a1 );
                        	                        t32 =-2.0 * ( p_surfQ[isurf][2] - a2 );
                                	                t42 =-2.0 * ( p_surfQ[isurf][3] - a3 );
                                        	        t52 =-2.0 * ( p_surfQ[isurf][4] - a4 );
						}
                                        }
					/*
				} else if (particle_on && link_point2[i][j][k][5]!=0 && link_point2[i][j][k][6]!=0) {
					isurf= p_surfi[i][j][k-1][0];
					jsurf= p_surfi[i][j][k  ][0];
					d2z0 = 4.0 * (p_surfQ[isurf][0] + p_surfQ[jsurf][0] - 2.0 * a0);
					d2z1 = 4.0 * (p_surfQ[isurf][1] + p_surfQ[jsurf][1] - 2.0 * a1);
					d2z2 = 4.0 * (p_surfQ[isurf][2] + p_surfQ[jsurf][2] - 2.0 * a2);
					d2z3 = 4.0 * (p_surfQ[isurf][3] + p_surfQ[jsurf][3] - 2.0 * a3);
					d2z4 = 4.0 * (p_surfQ[isurf][4] + p_surfQ[jsurf][4] - 2.0 * a4);
					t12  = p_surfQ[jsurf][0] - p_surfQ[isurf][0];
					t22  = p_surfQ[jsurf][1] - p_surfQ[isurf][1];
					t32  = p_surfQ[jsurf][2] - p_surfQ[isurf][2];
					t42  = p_surfQ[jsurf][3] - p_surfQ[isurf][3];
					t52  = p_surfQ[jsurf][4] - p_surfQ[isurf][4];
					 */
				} else {
					if(k==0 && wall_anchoring_on!=0){
						d2z0 = eight3rd*Qd[0] + four3rd*Qu[0] - 4.0 * a0;
						d2z1 = eight3rd*Qd[1] + four3rd*Qu[1] - 4.0 * a1;
						d2z2 = eight3rd*Qd[2] + four3rd*Qu[2] - 4.0 * a2;
						d2z3 = eight3rd*Qd[3] + four3rd*Qu[3] - 4.0 * a3;
						d2z4 = eight3rd*Qd[4] + four3rd*Qu[4] - 4.0 * a4;
						if(scheme==1) {
							t12 =-four3rd*Qd[0] + third*Qu[0] +   a0;
							t22 =-four3rd*Qd[1] + third*Qu[1] +   a1;
							t32 =-four3rd*Qd[2] + third*Qu[2] +   a2;
							t42 =-four3rd*Qd[3] + third*Qu[3] +   a3;
							t52 =-four3rd*Qd[4] + third*Qu[4] +   a4;
						} else if(scheme==2) {
							if(uz>0) {
								t12 = 2.0 * ( a0 - Qd[0] );
                	                                        t22 = 2.0 * ( a1 - Qd[1] );
                        	                                t32 = 2.0 * ( a2 - Qd[2] );
                                	                        t42 = 2.0 * ( a3 - Qd[3] );
                                        	                t52 = 2.0 * ( a4 - Qd[4] );
							} else if(uz<0) {
								t12 =-0.5 * Qu2[0] + 2.0 * Qu[0] - 1.5 * a0;
        	                                                t22 =-0.5 * Qu2[1] + 2.0 * Qu[1] - 1.5 * a1;
                	                                        t32 =-0.5 * Qu2[2] + 2.0 * Qu[2] - 1.5 * a2;
                        	                                t42 =-0.5 * Qu2[3] + 2.0 * Qu[3] - 1.5 * a3;
                                	                        t52 =-0.5 * Qu2[4] + 2.0 * Qu[4] - 1.5 * a4;
							}
						}
					} else if(k==Nz-1 && wall_anchoring_on!=0){
						d2z0 = eight3rd*Qu[0] + four3rd*Qd[0] - 4.0 * a0;
						d2z1 = eight3rd*Qu[1] + four3rd*Qd[1] - 4.0 * a1;
						d2z2 = eight3rd*Qu[2] + four3rd*Qd[2] - 4.0 * a2;
						d2z3 = eight3rd*Qu[3] + four3rd*Qd[3] - 4.0 * a3;
						d2z4 = eight3rd*Qu[4] + four3rd*Qd[4] - 4.0 * a4;
						if(scheme==1) {
							t12  = four3rd*Qu[0] - third*Qd[0] -  a0;
							t22  = four3rd*Qu[1] - third*Qd[1] -  a1;
							t32  = four3rd*Qu[2] - third*Qd[2] -  a2;
							t42  = four3rd*Qu[3] - third*Qd[3] -  a3;
							t52  = four3rd*Qu[4] - third*Qd[4] -  a4;
						} else if(scheme==2) {
	                                                if(uz<0) {
        	                                                t12 =-2.0 * ( a0 - Qu[0] );
                	                                        t22 =-2.0 * ( a1 - Qu[1] );
                        	                                t32 =-2.0 * ( a2 - Qu[2] );
                                	                        t42 =-2.0 * ( a3 - Qu[3] );
                                        	                t52 =-2.0 * ( a4 - Qu[4] );
	                                                } else if(uz<0) {
        	                                                t12 = 0.5 * Qd2[0] - 2.0 * Qd[0] + 1.5 * a0;
                	                                        t22 = 0.5 * Qd2[1] - 2.0 * Qd[1] + 1.5 * a1;
                        	                                t32 = 0.5 * Qd2[2] - 2.0 * Qd[2] + 1.5 * a2;
                                	                        t42 = 0.5 * Qd2[3] - 2.0 * Qd[3] + 1.5 * a3;
                                        	                t52 = 0.5 * Qd2[4] - 2.0 * Qd[4] + 1.5 * a4;
							}
                                                }
					} else {
						d2z0 = Qd[0] + Qu[0] - 2.0*a0;
						d2z1 = Qd[1] + Qu[1] - 2.0*a1;
						d2z2 = Qd[2] + Qu[2] - 2.0*a2;
						d2z3 = Qd[3] + Qu[3] - 2.0*a3;
						d2z4 = Qd[4] + Qu[4] - 2.0*a4;
						if(scheme==1) {
							t12  = 0.5 * ( Qu[0] - Qd[0] );
							t22  = 0.5 * ( Qu[1] - Qd[1] );
							t32  = 0.5 * ( Qu[2] - Qd[2] );
							t42  = 0.5 * ( Qu[3] - Qd[3] );
							t52  = 0.5 * ( Qu[4] - Qd[4] );
						} else if(scheme==2) {
							if(uz>0) {
								t12 = 1.5 * a0 - 2.0 * Qd[0] + 0.5 * Qd2[0];
                	                                        t22 = 1.5 * a1 - 2.0 * Qd[1] + 0.5 * Qd2[1];
                        	                                t32 = 1.5 * a2 - 2.0 * Qd[2] + 0.5 * Qd2[2];
                                	                        t42 = 1.5 * a3 - 2.0 * Qd[3] + 0.5 * Qd2[3];
                                        	                t52 = 1.5 * a4 - 2.0 * Qd[4] + 0.5 * Qd2[4];
							} else if(uz<0) {
        	                                                t12 =-1.5 * a0 + 2.0 * Qu[0] - 0.5 * Qu2[0];
                	                                        t22 =-1.5 * a1 + 2.0 * Qu[1] - 0.5 * Qu2[1];
                        	                                t32 =-1.5 * a2 + 2.0 * Qu[2] - 0.5 * Qu2[2];
                                	                        t42 =-1.5 * a3 + 2.0 * Qu[3] - 0.5 * Qu2[3];
                                        	                t52 =-1.5 * a4 + 2.0 * Qu[4] - 0.5 * Qu2[4];
							}
						}
					}
				}
				
				d2Q[0][0] = d2x0 + d2y0 + d2z0;
				d2Q[0][1] = d2x1 + d2y1 + d2z1;
				d2Q[0][2] = d2x2 + d2y2 + d2z2;
				d2Q[1][1] = d2x3 + d2y3 + d2z3;
				d2Q[1][2] = d2x4 + d2y4 + d2z4;
				d2Q[1][0] = d2Q[0][1];
				d2Q[2][0] = d2Q[0][2];
				d2Q[2][1] = d2Q[1][2];
				d2Q[2][2] =-d2Q[0][0] - d2Q[1][1];

//				lap_Q[i][j][k]=2.0*(d2Q[0][1]+d2Q[0][2]+d2Q[1][2]);
				
//	calculate H				
				trQ2 = 2.0*(a0*a0+a1*a1+a2*a2+a3*a3+a4*a4+a0*a3);

				for(ii=0;ii<3;ii++){
					for(jj=0;jj<3;jj++){
						temp = Q0[i][j][k][ii][jj];
						H0[i][j][k][ii][jj] = -1*(A_ldg*(1.0-third*U)*temp - A_ldg*U*( QQ(Q0[i][j][k],ii,jj) - trQ2*(temp + OneThirdDelta(ii,jj))) - kappa * d2Q[ii][jj]);
					}
				}

				if(debug_on!=0){
					trH3rd = third*(H0[i][j][k][0][0]+H0[i][j][k][1][1]+H0[i][j][k][2][2]);
					if (trH3rd>1e-15 || trH3rd<-1e-15)printf("trH3rd nonzero: %d %d %d: %le\n",i, j, k,trH3rd);
				}
				
//	calculate convective Q
				ux = u[i][j][k][0];
				uy = u[i][j][k][1];
				uz = u[i][j][k][2];
				convQ[i][j][k][0] = ux * t10 + uy * t11 + uz * t12;
				convQ[i][j][k][1] = ux * t20 + uy * t21 + uz * t22;
				convQ[i][j][k][2] = ux * t30 + uy * t31 + uz * t32;
				convQ[i][j][k][3] = ux * t40 + uy * t41 + uz * t42;
				convQ[i][j][k][4] = ux * t50 + uy * t51 + uz * t52;

			}
		}
	}
	
}




// calculate Q tensor related sigma
void cal_stress_q(real Q0[Nx][Ny][Nz][3][3], real H0[Nx][Ny][Nz][3][3], real sigma_q[Nx][Ny][Nz][3][3])
{

	int i, j, k, ii=0, jj=0, kk=0, mm;
	int ip, im, jp, jm, kp, km;
	int ipar, isurf, jsurf;
	real t10, t11, t12, t20, t21, t22, t30, t31, t32, t40, t41, t42, t50, t51, t52;
	real temp, Qu[5], Qd[5];
	real dx, dy, dz, dr2i;
	real a0, a1, a2, a3, a4;
	real trQ2, trH3rd, d2Q[3][3], dQ2, dQ[3][3];
	real sum;
	
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			for (k=0; k<Nz; k++) {
				a0 = Q0[i][j][k][0][0];
				a1 = Q0[i][j][k][0][1];
				a2 = Q0[i][j][k][0][2];
				a3 = Q0[i][j][k][1][1];
				a4 = Q0[i][j][k][1][2];
				
		    	im = i-1;
		  		jm = j-1;
				ip = i+1;
		  	    jp = j+1;
				km = k-1;
				kp = k+1;

				if(i==0) im = Nx-1;
		  		if(j==0) jm = Ny-1;
		  		if(k==0) km = Nz-1;

	     		if(i==Nx-1) ip = 0;
		  		if(j==Ny-1) jp = 0;
		  		if(k==Nz-1) kp = 0;
			
				if(k==0 && wall_anchoring_on!=0){
					Qd[0] = Qsurf[i][j][0][0];
					Qd[1] = Qsurf[i][j][0][1];
					Qd[2] = Qsurf[i][j][0][2];
					Qd[3] = Qsurf[i][j][0][3];
					Qd[4] = Qsurf[i][j][0][4];
				} else {
					Qd[0] = Q0[i][j][km][0][0];
					Qd[1] = Q0[i][j][km][0][1];
					Qd[2] = Q0[i][j][km][0][2];
					Qd[3] = Q0[i][j][km][1][1];
					Qd[4] = Q0[i][j][km][1][2];
				}
				if(k==Nz-1 && wall_anchoring_on!=0) {
					Qu[0] = Qsurf[i][j][1][0];
					Qu[1] = Qsurf[i][j][1][1];
					Qu[2] = Qsurf[i][j][1][2];
					Qu[3] = Qsurf[i][j][1][3];
					Qu[4] = Qsurf[i][j][1][4];
				} else {
					Qu[0] = Q0[i][j][kp][0][0];
					Qu[1] = Q0[i][j][kp][0][1];
					Qu[2] = Q0[i][j][kp][0][2];
					Qu[3] = Q0[i][j][kp][1][1];
					Qu[4] = Q0[i][j][kp][1][2];
				}
				
				ipar=0;
				if(particle_on && link_point2[i][j][k][1]!=0 && link_point2[i][j][k][3]==0) {
					isurf= p_surfi[i][j][k][0];
					t10  = four3rd*p_surfQ[isurf][0]  - third*Q0[im][j][k][0][0]   -      a0;
					t20  = four3rd*p_surfQ[isurf][1]  - third*Q0[im][j][k][0][1]   -      a1;
					t30  = four3rd*p_surfQ[isurf][2]  - third*Q0[im][j][k][0][2]   -      a2;
					t40  = four3rd*p_surfQ[isurf][3]  - third*Q0[im][j][k][1][1]   -      a3;
					t50  = four3rd*p_surfQ[isurf][4]  - third*Q0[im][j][k][1][2]   -      a4;
				} else if(particle_on && link_point2[i][j][k][3]!=0 && link_point2[i][j][k][1]==0) {
					isurf= p_surfi[i-1][j][k][0];
					t10  =-four3rd*p_surfQ[isurf][0]  + third*Q0[ip][j][k][0][0]   +      a0;
					t20  =-four3rd*p_surfQ[isurf][1]  + third*Q0[ip][j][k][0][1]   +      a1;
					t30  =-four3rd*p_surfQ[isurf][2]  + third*Q0[ip][j][k][0][2]   +      a2;
					t40  =-four3rd*p_surfQ[isurf][3]  + third*Q0[ip][j][k][1][1]   +      a3;
					t50  =-four3rd*p_surfQ[isurf][4]  + third*Q0[ip][j][k][1][2]   +      a4;
					/*
				} else if (particle_on && link_point2[i][j][k][1]!=0 && link_point2[i][j][k][3]!=0) {
					isurf= p_surfi[i-1][j][k][0];
					jsurf= p_surfi[i  ][j][k][0];
					t10  = p_surfQ[jsurf][0] - p_surfQ[isurf][0];
					t20  = p_surfQ[jsurf][1] - p_surfQ[isurf][1];
					t30  = p_surfQ[jsurf][2] - p_surfQ[isurf][2];
					t40  = p_surfQ[jsurf][3] - p_surfQ[isurf][3];
					t50  = p_surfQ[jsurf][4] - p_surfQ[isurf][4];
					 */
				} else {
					t10  = 0.5 * (Q0[ip][j][k][0][0]-Q0[im][j][k][0][0]);
					t20  = 0.5 * (Q0[ip][j][k][0][1]-Q0[im][j][k][0][1]);
					t30  = 0.5 * (Q0[ip][j][k][0][2]-Q0[im][j][k][0][2]);
					t40  = 0.5 * (Q0[ip][j][k][1][1]-Q0[im][j][k][1][1]);
					t50  = 0.5 * (Q0[ip][j][k][1][2]-Q0[im][j][k][1][2]);
				}
				
				if(particle_on && link_point2[i][j][k][2]!=0 && link_point2[i][j][k][4]==0) {
					isurf= p_surfi[i][j][k][1];
					t11  = four3rd*p_surfQ[isurf][0]  - third*Q0[i][jm][k][0][0]   -      a0;
					t21  = four3rd*p_surfQ[isurf][1]  - third*Q0[i][jm][k][0][1]   -      a1;
					t31  = four3rd*p_surfQ[isurf][2]  - third*Q0[i][jm][k][0][2]   -      a2;
					t41  = four3rd*p_surfQ[isurf][3]  - third*Q0[i][jm][k][1][1]   -      a3;
					t51  = four3rd*p_surfQ[isurf][4]  - third*Q0[i][jm][k][1][2]   -      a4;
				} else if(particle_on && link_point2[i][j][k][4]!=0 && link_point2[i][j][k][2]==0) {
					isurf= p_surfi[i][j-1][k][1];
					t11  =-four3rd*p_surfQ[isurf][0]  + third*Q0[i][jp][k][0][0]   +      a0;
					t21  =-four3rd*p_surfQ[isurf][1]  + third*Q0[i][jp][k][0][1]   +      a1;
					t31  =-four3rd*p_surfQ[isurf][2]  + third*Q0[i][jp][k][0][2]   +      a2;
					t41  =-four3rd*p_surfQ[isurf][3]  + third*Q0[i][jp][k][1][1]   +      a3;
					t51  =-four3rd*p_surfQ[isurf][4]  + third*Q0[i][jp][k][1][2]   +      a4;
					/*
				} else if (particle_on && link_point2[i][j][k][2]!=0 && link_point2[i][j][k][4]!=0) {
					isurf= p_surfi[i][j-1][k][0];
					jsurf= p_surfi[i][j  ][k][0];
					t11  = p_surfQ[jsurf][0] - p_surfQ[isurf][0];
					t21  = p_surfQ[jsurf][1] - p_surfQ[isurf][1];
					t31  = p_surfQ[jsurf][2] - p_surfQ[isurf][2];
					t41  = p_surfQ[jsurf][3] - p_surfQ[isurf][3];
					t51  = p_surfQ[jsurf][4] - p_surfQ[isurf][4];
					 */
				} else {
					t11  = 0.5 * (Q0[i][jp][k][0][0]-Q0[i][jm][k][0][0]);
					t21  = 0.5 * (Q0[i][jp][k][0][1]-Q0[i][jm][k][0][1]);
					t31  = 0.5 * (Q0[i][jp][k][0][2]-Q0[i][jm][k][0][2]);
					t41  = 0.5 * (Q0[i][jp][k][1][1]-Q0[i][jm][k][1][1]);
					t51  = 0.5 * (Q0[i][jp][k][1][2]-Q0[i][jm][k][1][2]);
				}
				
				if(particle_on && link_point2[i][j][k][5]!=0 && link_point2[i][j][k][6]==0) {
					isurf= p_surfi[i][j][k][2];
					t12  = four3rd*p_surfQ[isurf][0]  - third * Qd[0]   -       a0;
					t22  = four3rd*p_surfQ[isurf][1]  - third * Qd[1]   -       a1;
					t32  = four3rd*p_surfQ[isurf][2]  - third * Qd[2]   -       a2;
					t42  = four3rd*p_surfQ[isurf][3]  - third * Qd[3]   -       a3;
					t52  = four3rd*p_surfQ[isurf][4]  - third * Qd[4]   -       a4;
				} else if(particle_on && link_point2[i][j][k][6]!=0 && link_point2[i][j][k][5]==0) {
					isurf= p_surfi[i][j][k-1][2];
					t12  =-four3rd*p_surfQ[isurf][0]  + third * Qu[0]   +       a0;
					t22  =-four3rd*p_surfQ[isurf][1]  + third * Qu[1]   +       a1;
					t32  =-four3rd*p_surfQ[isurf][2]  + third * Qu[2]   +       a2;
					t42  =-four3rd*p_surfQ[isurf][3]  + third * Qu[3]   +       a3;
					t52  =-four3rd*p_surfQ[isurf][4]  + third * Qu[4]   +       a4;
					/*
				} else if (particle_on && link_point2[i][j][k][5]!=0 && link_point2[i][j][k][6]!=0) {
					isurf= p_surfi[i][j][k-1][0];
					jsurf= p_surfi[i][j][k  ][0];
					t12  = p_surfQ[jsurf][0] - p_surfQ[isurf][0];
					t22  = p_surfQ[jsurf][1] - p_surfQ[isurf][1];
					t32  = p_surfQ[jsurf][2] - p_surfQ[isurf][2];
					t42  = p_surfQ[jsurf][3] - p_surfQ[isurf][3];
					t52  = p_surfQ[jsurf][4] - p_surfQ[isurf][4];
					 */
				} else {
					if(k==0 && wall_anchoring_on!=0){
						t12  =-four3rd*Qd[0] + third*Qu[0]    +       a0;
						t22  =-four3rd*Qd[1] + third*Qu[1]    +       a1;
						t32  =-four3rd*Qd[2] + third*Qu[2]    +       a2;
						t42  =-four3rd*Qd[3] + third*Qu[3]    +       a3;
						t52  =-four3rd*Qd[4] + third*Qu[4]    +       a4;
					} else if(k==Nz-1 && wall_anchoring_on!=0){
						t12  = four3rd*Qu[0] - third*Qd[0]    -       a0;
						t22  = four3rd*Qu[1] - third*Qd[1]    -       a1;
						t32  = four3rd*Qu[2] - third*Qd[2]    -       a2;
						t42  = four3rd*Qu[3] - third*Qd[3]    -       a3;
						t52  = four3rd*Qu[4] - third*Qd[4]    -       a4;
					} else {
						t12  = 0.5 * ( Qu[0] - Qd[0] );
						t22  = 0.5 * ( Qu[1] - Qd[1] );
						t32  = 0.5 * ( Qu[2] - Qd[2] );
						t42  = 0.5 * ( Qu[3] - Qd[3] );
						t52  = 0.5 * ( Qu[4] - Qd[4] );
					}
				}	
				
				dQ[0][0]=2.*(t10*t10+t20*t20+t30*t30+t40*t40+t50*t50+t10*t40);
				dQ[1][1]=2.*(t11*t11+t21*t21+t31*t31+t41*t41+t51*t51+t11*t41);
				dQ[2][2]=2.*(t12*t12+t22*t22+t32*t32+t42*t42+t52*t52+t12*t42);
				dQ[0][1]=2.*(t10*t11+t20*t21+t30*t31+t40*t41+t50*t51)+t10*t41+t11*t40;
				dQ[0][2]=2.*(t10*t12+t20*t22+t30*t32+t40*t42+t50*t52)+t10*t42+t12*t40;
				dQ[1][2]=2.*(t11*t12+t21*t22+t31*t32+t41*t42+t51*t52)+t11*t42+t12*t41;
				dQ[1][0]=dQ[0][1];
				dQ[2][0]=dQ[0][2];
				dQ[2][1]=dQ[1][2];
				
				dQ2 = dQ[0][0]+dQ[1][1]+dQ[2][2];
//				lap_Q[i][j][k] = dQ2;

				for(sum=0.0,ii=0;ii<3;ii++){
					for(jj=0;jj<3;jj++){
						sum = sum + Q0[i][j][k][ii][jj]*H0[i][j][k][ii][jj];
					}
		  		}

		  		for(ii=0;ii<3;ii++){
					for(jj=0;jj<3;jj++){
						temp = 0.5*kappa*dQ2*Delta(ii,jj) + 2*xi*(Q0[i][j][k][ii][jj] + OneThirdDelta(ii,jj))*sum;
//						temp = 2.0*xi*(Q0[i][j][k][ii][jj] + OneThirdDelta(ii,jj))*sum;
						for(mm=0;mm<3;mm++){
							temp += - xi*H0[i][j][k][ii][mm]*(Q0[i][j][k][mm][jj] + OneThirdDelta(mm,jj)) - xi*(Q0[i][j][k][ii][mm] + OneThirdDelta(ii,mm))*H0[i][j][k][mm][jj] ;
						}
						for(mm=0;mm<3;mm++){
							temp+=Q0[i][j][k][ii][mm]*H0[i][j][k][mm][jj]-H0[i][j][k][ii][mm]*Q0[i][j][k][mm][jj];
						}

						temp += -kappa * dQ[ii][jj];
						if (t_current>-1){
							sigma_q[i][j][k][ii][jj] = temp;
							if (ii==jj ) {
//								sigma_q[i][j][k][ii][jj] = 0;
							}
						} else {
							sigma_q[i][j][k][ii][jj] = 0;	// decouple velocity
						}
					}
				}

				
				
			}
		}
	}
	
//	printf("\nsigma_q=%e %e %e %e %e %e %e %e %e\n",sigma_q[0][0][Nz/2][0][0],sigma_q[0][0][Nz/2][0][1],sigma_q[0][0][Nz/2][0][2],sigma_q[0][0][Nz/2][1][0],sigma_q[0][0][Nz/2][1][1],sigma_q[0][0][Nz/2][1][2],sigma_q[0][0][Nz/2][2][0],sigma_q[0][0][Nz/2][2][1],sigma_q[0][0][Nz/2][2][2]);
}




void cal_sigma_p(real Q0[Nx][Ny][Nz][3][3])
{

	int i, j, k, ii=0, jj=0, kk=0;
	int ip, im, jp, jm, kp, km;
	int ipar, isurf, jsurf;
	real t10, t11, t12, t20, t21, t22, t30, t31, t32, t40, t41, t42, t50, t51, t52;
	real temp, Qu[5], Qd[5];
	real d2x0, d2x1, d2x2, d2x3, d2x4, d2y0, d2y1, d2y2, d2y3, d2y4, d2z0, d2z1, d2z2, d2z3, d2z4;
	real dx, dy, dz, dr2i;
	real a0, a1, a2, a3, a4;
	real trQ2, trH3rd, d2Q[3][3], dQ2;
	real d2Q1, d2Q2, d2Q3, d2Q4, d2Q5;
	real ux, uy, uz;

/*
	for (i=0; i<Nx; i++) {
		for (j=0; j<Ny; j++) {
			sigma_p[i][j][0][0] = 0;
			sigma_p[i][j][0][2] = 0;
			sigma_p[i][j][0][1] = sigma_q[i][j][1][1][2] - sigma_q[i][j][0][1][2];
			for (k=1; k<Nz-1; k++) {
				sigma_p[i][j][k][0] = 0;
				sigma_p[i][j][k][2] = 0;
				sigma_p[i][j][k][1] = 0.5 * (sigma_q[i][j][k+1][1][2] - sigma_q[i][j][k-1][1][2]);
			}
			sigma_p[i][j][Nz-1][0] = 0;
			sigma_p[i][j][Nz-1][2] = 0;
			sigma_p[i][j][Nz-1][1] = sigma_q[i][j][Nz-1][1][2] - sigma_q[i][j][Nz-2][1][2];
		}
	}
*/
	

	for (i=0; i<Nx; i++) {
		im = i - 1;
		ip = i + 1;
		if(i==0) im = Nx-1;
		if(i==Nx-1) ip = 0;
		for (j=0; j<Ny; j++) {
			jm = j - 1;
			jp = j + 1;
			if(j==0) jm = Ny-1;
			if(j==Ny-1) jp = 0;
			for (k=0; k<Nz; k++) {
				km = k - 1;
				kp = k + 1;
				if(k==0) km = Nz-1;
				if(k==Nz-1) kp = 0;
				sigma_p[i][j][k][0] = 0.5 * ( sigma_q[ip][j][k][0][0] - sigma_q[im][j][k][0][0] + sigma_q[i][jp][k][0][1] - sigma_q[i][jm][k][0][1] );
				sigma_p[i][j][k][1] = 0.5 * ( sigma_q[ip][j][k][1][0] - sigma_q[im][j][k][1][0] + sigma_q[i][jp][k][1][1] - sigma_q[i][jm][k][1][1] );
				sigma_p[i][j][k][2] = 0.5 * ( sigma_q[ip][j][k][2][0] - sigma_q[im][j][k][2][0] + sigma_q[i][jp][k][2][1] - sigma_q[i][jm][k][2][1] );
				if ( k==0 && wall_on!=0) {
					sigma_p[i][j][k][0] += -0.5*sigma_q[i][j][2][0][2] + 2.0*sigma_q[i][j][1][0][2] - 1.5*sigma_q[i][j][0][0][2];
//					sigma_p[i][j][k][0] += sigma_q[i][j][1][0][2] - sigma_q[i][j][0][0][2];
					sigma_p[i][j][k][1] += -0.5*sigma_q[i][j][2][1][2] + 2.0*sigma_q[i][j][1][1][2] - 1.5*sigma_q[i][j][0][1][2];
					sigma_p[i][j][k][2] += -0.5*sigma_q[i][j][2][2][2] + 2.0*sigma_q[i][j][1][2][2] - 1.5*sigma_q[i][j][0][2][2];
				} else if ( k==Nz-1 && wall_on!=0) {
					sigma_p[i][j][k][0] += 0.5*sigma_q[i][j][Nz-3][0][2] - 2.0*sigma_q[i][j][Nz-2][0][2] + 1.5*sigma_q[i][j][Nz-1][0][2];
//					sigma_p[i][j][k][0] += sigma_q[i][j][Nz-1][0][2] - sigma_q[i][j][Nz-2][0][2];
					sigma_p[i][j][k][1] += 0.5*sigma_q[i][j][Nz-3][1][2] - 2.0*sigma_q[i][j][Nz-2][1][2] + 1.5*sigma_q[i][j][Nz-1][1][2];
					sigma_p[i][j][k][2] += 0.5*sigma_q[i][j][Nz-3][2][2] - 2.0*sigma_q[i][j][Nz-2][2][2] + 1.5*sigma_q[i][j][Nz-1][2][2];
				} else {
					sigma_p[i][j][k][0] += 0.5 * ( sigma_q[i][j][kp][0][2] - sigma_q[i][j][km][0][2] );
					sigma_p[i][j][k][1] += 0.5 * ( sigma_q[i][j][kp][1][2] - sigma_q[i][j][km][1][2] );
					sigma_p[i][j][k][2] += 0.5 * ( sigma_q[i][j][kp][2][2] - sigma_q[i][j][km][2][2] );
				}
			}
		}
	}

 
}
