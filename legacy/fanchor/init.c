#include "lb.h"

void init(real f[Nx][Ny][Nz][15], real Q0[Nx][Ny][Nz][3][3])
{
	int i,j,k,ii,  ip, jp, kp, im, jm, km;
	double sita, sita0, lambda1, lambda2, nlocal[3], slocal, slocali;

	sita0 = ntop[0]*nbot[0]+ntop[1]*nbot[1]+ntop[2]*nbot[2];
	sita0 = sita0/sqrt(ntop[0]*ntop[0]+ntop[1]*ntop[1]+ntop[2]*ntop[2]);        
	sita0 = sita0/sqrt(nbot[0]*nbot[0]+nbot[1]*nbot[1]+nbot[2]*nbot[2]);
	sita0 = acos(sita0);

	if(newrun_on!=0){
		for (i=0; i<Nx; i++) {
			for (j=0; j<Ny; j++) {
				for (k=0; k<Nz; k++) {
					Rho[i][j][k]=rho;
					u[i][j][k][0]=0;
					u[i][j][k][1]=uy_bottom+(uy_top-uy_bottom)*((double)k+0.5)/((double)Nz);
					u[i][j][k][2]=0;
					
					clear2d(Q0[i][j][k]);

					if (sita0<1e-2 && sita0>-1e-2) {
						lambda1 = 0.5;
						lambda2 = 0.5;
					} else {
						sita = sita0 * ((double)k+0.5)/((double)Nz);
						lambda1=(cos(sita)-cos(sita0)*cos(sita0-sita));
						lambda2=(cos(sita0-sita)-cos(sita)*cos(sita0));
					}
					nlocal[0] = lambda1 * nbot[0] + lambda2 * ntop[0];
					nlocal[1] = lambda1 * nbot[1] + lambda2 * ntop[1];
					nlocal[2] = lambda1 * nbot[2] + lambda2 * ntop[2];
					slocal = nlocal[0]*nlocal[0] + nlocal[1]*nlocal[1] + nlocal[2]*nlocal[2];
					slocali= 1.0/slocal;
					Q0[i][j][k][0][0] =  S*(nlocal[0]*nlocal[0]*slocali - third);
					Q0[i][j][k][0][1] =  S*(nlocal[0]*nlocal[1]*slocali);
					Q0[i][j][k][0][2] =  S*(nlocal[0]*nlocal[2]*slocali);
					Q0[i][j][k][1][1] =  S*(nlocal[1]*nlocal[1]*slocali - third);
					Q0[i][j][k][1][2] =  S*(nlocal[1]*nlocal[2]*slocali);
					Q0[i][j][k][1][0] =  Q0[i][j][k][0][1];
					Q0[i][j][k][2][0] =  Q0[i][j][k][0][2];
					Q0[i][j][k][2][1] =  Q0[i][j][k][1][2];
					Q0[i][j][k][2][2] = -Q0[i][j][k][0][0]-Q0[i][j][k][1][1];
				}
			}
		}
	}
	else {
		read_restart(Q0,Rho,u);
	}
//	printf("rho=%25.20f, uz=%25.20f",Rho[0][0][0],u[0][0][0][2]);

	                //For finite anchoring
                for(i=0;i<Nx;i++){
                        for(j=0;j<Ny;j++){

                                Qsurf[i][j][0][0] =  S*(nbot[0]*nbot[0] - third);
                                Qsurf[i][j][0][1] =  S*(nbot[0]*nbot[1]);
                                Qsurf[i][j][0][2] =  S*(nbot[0]*nbot[2]);
                                Qsurf[i][j][0][3] =  S*(nbot[1]*nbot[1] - third);
                                Qsurf[i][j][0][4] =  S*(nbot[1]*nbot[2]);

                                Qsurf[i][j][1][0] =  S*(ntop[0]*ntop[0] - third);
                                Qsurf[i][j][1][1] =  S*(ntop[0]*ntop[1]);
                                Qsurf[i][j][1][2] =  S*(ntop[0]*ntop[2]);
                                Qsurf[i][j][1][3] =  S*(ntop[1]*ntop[1] - third);
                                Qsurf[i][j][1][4] =  S*(ntop[1]*ntop[2]);

                                Qsurf0[i][j][0][0]=  S*(nbot[0]*nbot[0] - third);
                                Qsurf0[i][j][0][1]=  S*(nbot[0]*nbot[1]);
                                Qsurf0[i][j][0][2]=  S*(nbot[0]*nbot[2]);
                                Qsurf0[i][j][0][3]=  S*(nbot[1]*nbot[1] - third);
                                Qsurf0[i][j][0][4]=  S*(nbot[1]*nbot[2]);

                                Qsurf0[i][j][1][0]=  S*(ntop[0]*ntop[0] - third);
                                Qsurf0[i][j][1][1]=  S*(ntop[0]*ntop[1]);
                                Qsurf0[i][j][1][2]=  S*(ntop[0]*ntop[2]);
                                Qsurf0[i][j][1][3]=  S*(ntop[1]*ntop[1] - third);
                                Qsurf0[i][j][1][4]=  S*(ntop[1]*ntop[2]);
//				for(k=0;k<Nz;k++){
//					u[i][j][k][0]=0;
//                                      u[i][j][k][1]=uy_bottom+(uy_top-uy_bottom)*(double)k/((double)(Nz-1));
//                                      u[i][j][k][2]=0;
//				}
                        }
                }
	

	cal_dQ(Q0, H, convQ);
	if(newrun_on==0)cal_stress_q(Q0, H, sigma_q);
	cal_sigma(Q0);
	cal_fequ(f,Rho,u,sigma);
	cal_W(u);

	if(newrun_on!=0){
		for(i=0;i<Nx;i++){
		    for(j=0;j<Ny;j++){
		    	for(k=0;k<Nz;k++){
		    		sigma_p[i][j][k][0] = 0;
		    		sigma_p[i][j][k][1] = 0;
		    		sigma_p[i][j][k][2] = 0;
		    	}
		    }
		}
	} else{
		cal_sigma_p(Q0);
	}
}
