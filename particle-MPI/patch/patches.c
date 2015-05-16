#include "main.h"
#include "particle.h"

void add_patch()
{
	FILE *param;
	int i, j, id, imin, imax, jmin, jmax, ii;
	double x, y, z, q[6];

	if (patch_on!=0) {
                param = fopen("patch.in", "r");

                fscanf(param, "patch_length %lf %lf %lf\n", &patchx, &patchy, &patchz);
                fscanf(param, "type W_patch %d %lf\n", &type_patch, &W_patch);
                fscanf(param, "n_patch %lf %lf %lf\n", &x, &y, &z);
                fclose(param);
	}

	if (patch_on!=0 && wall_on!=0 && Q_on!=0 && myid==0) {
		imin = ceil(((double)Nx - 1.0 - patchx) * 0.5);
		imax =      ((double)Nx - 1.0 + patchx) * 0.5;
		jmin = ceil(((double)Ny - 1.0 - patchy) * 0.5);
		jmax =      ((double)Ny - 1.0 + patchy) * 0.5;

                if (imin<0)   imin = 0;
                if (imax>=Nx) imax =Nx-1;
		if (jmin<0)   jmin = 0;
		if (jmax>=Ny) jmax =Ny-1;

		for (j=jmin; j<=jmax; j++) {
			for (i=imin; i<=imax; i++) {
				id = i + j*Nx;
				surf[id*10]  = type_patch;
				surf[id*10+1]= W_patch;
				ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
				for (ii=0; ii<5; ii++) surf[id*10+5+ii]=q[ii];
			}
		}
	}
}

int onpatch(int idx, int idy, int idz)
{
	int imin, imax, jmin, jmax, kmax;

        if (patch_on!=0) {
		imin = ceil(((double)Nx - 1.0 - patchx) * 0.5);
		imax =      ((double)Nx - 1.0 + patchx) * 0.5;
		jmin = ceil(((double)Ny - 1.0 - patchy) * 0.5);
		jmax =      ((double)Ny - 1.0 + patchy) * 0.5;
		kmax = patchz;
		if (idx>=imin && idx<=imax && idy>=jmin && idy <=jmax && idz<=kmax) return 1; 
	}
	return -1;
}

void patch_iden()
{
	int xlo, xhi, ylo, yhi, zlo, zhi, x1, y1, z1, x2, y2, z2, i1, j1, k1, i2, j2, k2, i3, j3, k3;
	int i, j, k, ii, id, id1, id2, id3, di, p1, p2, nlink=0;
	int si, sj, sk, sid, sii, sjj;
	double r, x, y, z, ub_dot_e=0., q[6];
	
	//		boundaries of search box 
	xlo = ceil(((double)Nx - 1.0 - patchx) * 0.5)-1;
	xhi =      ((double)Nx - 1.0 + patchx) * 0.5 +1;
	ylo = ceil(((double)Ny - 1.0 - patchy) * 0.5)-1;
	yhi =      ((double)Ny - 1.0 + patchy) * 0.5 +1;
	zlo = 0;
	zhi = patchz+1;
	if (wall_on!=0) {
		if(zlo <  0)  zlo = 0;
		if(zhi >= Nz) zhi = Nz-1;
	}
	for (z1=zlo; z1<=zhi; z1++) {
		for (y1=ylo; y1<=yhi; y1++) {
			for (x1=xlo; x1<=xhi; x1++) {
				p1=onpatch(x1, y1, z1);
				i1 = x1;
				if(x1<0)   i1 += Nx;
				if(x1>=Nx) i1 -= Nx;
				j1 = y1;
				if(y1<0)   j1 += Ny;
				if(y1>=Ny) j1 -= Ny;
				k1 = z1;
				if (wall_on==0) {
					if(k1<0)   k1 += Nz;
					if(k1>=Nz) k1 -= Nz;
				}
				id1 = i1 + (j1+k1*Ny)*Nx;
				//	inside patch
				if (p1>0 && id1/point==myid) info[id1%point]=-2;
				
				for (j=2; j<=10; j++) {
					z2 = z1 + e[j][2];
					if ( j!=3 && j!=5 && ( wall_on==0 || (z2>=0 && z2<Nz) ) ) {
						x2 = x1 + e[j][0];
						y2 = y1 + e[j][1];
						p2=onpatch(x2, y2, z2);
						
						if (p1 * p2 < 0) {
							i2 = x2;
							if(i2<0)   i2 += Nx;
							if(i2>=Nx) i2 -= Nx;
							j2 = y2;
							if(j2<0)   j2 += Ny;
							if(j2>=Ny) j2 -= Ny;
							k2 = z2;
							if (wall_on==0) {
								if(k2<0)   k2 += Nz;
								if(k2>=Nz) k2 -= Nz;
							}
							k   = bounce(j);
							id2 = i2 + (j2+k2*Ny)*Nx;
							
							x   = (double)x1 - (double)x2;
							y   = (double)y1 - (double)y2;
							z   = (double)z1 - (double)z2;
							if (p1>0) {
								x = -x;
								y = -y;
								z = -z;
							}
							
							if (flow_on!=0) {			//	modify streaming
								ub_dot_e = ux_bot*e[j][0] + uy_bot*e[j][1];
								if(id1/point==myid) {
									nextf[(id1%point)*15+j]=(id1%point)*15+k;
								}
								if(id2/point==myid) {
									nextf[(id2%point)*15+k]=(id2%point)*15+j;
								}
								nlink++;
							}
							
							if (j<7 && Q_on!=0) {			//	create surface point
								if (nsurf/node==myid) {
									id        = 10*(nsurf%node);
									surf[id]  = type_patch;
									surf[id+1]= W_patch;
									normalize(&x,&y,&z);
									surf[id+2]= x;
									surf[id+3]= y;
									surf[id+4]= z;
									if (type_patch<=1) {
	                                                                        ntoq(x, y, z, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
									} else {
	                                                                        ntoq(y, z, x, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
									}
									for (ii=0; ii<5; ii++) {
										surf[id+5+ii]  = q[ii];
									}
								}
								
								//	look for the nearest bulk point for the surface point
								if (p1<0) {
									si = i1;
									sj = j1;
									sk = k1;
									sid= id1;
									sii= j;
									sjj= k;
								} else {
									si = i2;
									sj = j2;
									sk = k2;
									sid= id2;
									sii= k;
									sjj= j;
								}
								if(sid/point==myid) neighb[(sid%point)*6+sii-1]=(nsurf-myid*node)*5-1;
								
								//	do the first pair neighbors
								i3 = si + e[sjj][0];
								j3 = sj + e[sjj][1];
								k3 = sk + e[sjj][2];
								if(i3<0)   i3 += Nx;
								if(i3>=Nx) i3 -= Nx;
								if(j3<0)   j3 += Ny;
								if(j3>=Ny) j3 -= Ny;
								if (wall_on==0) {
									if(k3<0)   k3 += Nz;
									if(k3>=Nz) k3 -= Nz;
								}
								id3 = i3 +(j3+k3*Ny)*Nx;
								if(nsurf/node==myid) {
									if (sii<sjj) {
										neighbsurf[(nsurf%node)*6+sii-1]=(sid-myid*point)*5;
										neighbsurf[(nsurf%node)*6+sjj-1]=(id3-myid*point)*5;
									} else {
										neighbsurf[(nsurf%node)*6+sjj-1]=(sid-myid*point)*5;
										neighbsurf[(nsurf%node)*6+sii-1]=(id3-myid*point)*5;
									}
									
									// find neighbors for surface point
									if (sii!=1 && sjj!=1) {
										if (x>0) di= 1;
										if (x<0) di=-1;
										i3 = si - di;
										j3 = sj;
										k3 = sk;
										if(i3<0)   i3 += Nx;
										if(i3>=Nx) i3 -= Nx;
										
										if (onpatch(i3,j3,k3)==-1) {
											neighbsurf[(nsurf%node)*6]  = 5*(i3 +(j3+k3*Ny)*Nx - myid*point)-1;
											i3 = si + di;
											if(i3<0)   i3 += Nx;
											if(i3>=Nx) i3 -= Nx;
											neighbsurf[(nsurf%node)*6+1]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point);
										} else {
											i3 = si + di;
											if(i3<0)   i3 += Nx;
											if(i3>=Nx) i3 -= Nx;
											neighbsurf[(nsurf%node)*6]=   5*(i3 +(j3+k3*Ny)*Nx - myid*point)-2;
											i3 = si + 2*di;
											if(i3<0)   i3 += Nx;
											if(i3>=Nx) i3 -= Nx;
											neighbsurf[(nsurf%node)*6+1]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point);
										}
									}
									
									if (sii!=3 && sjj!=3) {
										if (y>0) di= 1;
										if (y<0) di=-1;
										i3 = si;
										j3 = sj - di;
										k3 = sk;
										if(j3<0)   j3 += Ny;
										if(j3>=Ny) j3 -= Ny;
										
										if (onpatch(i3,j3,k3)==-1) {
											neighbsurf[(nsurf%node)*6+2]= 5*(i3 +(j3+k3*Ny)*Nx -myid*point)-1;
											j3 = sj + di;
											if(j3<0)   j3 += Ny;
											if(j3>=Nx) j3 -= Ny;
											neighbsurf[(nsurf%node)*6+3]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point);
										} else {
											j3 = sj + di;
											if(j3<0)   j3 += Ny;
											if(j3>=Nx) j3 -= Ny;
											neighbsurf[(nsurf%node)*6+2]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point)-2;
											j3 = sj + 2*di;
											if(j3<0)   j3 += Ny;
											if(j3>=Nx) j3 -= Ny;
											neighbsurf[(nsurf%node)*6+3]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point);
										}
									}
									
									if (sii!=5 && sjj!=5) {
										if (z>0) di= 1;
										if (z<0) di=-1;
										i3 = si;
										j3 = sj;
										k3 = sk - di;
										if (wall_on==0) {
											if(k3<0)   k3 += Nz;
											if(k3>=Nz) k3 -= Nz;
										}
										
										if (onpatch(i3,j3,k3)==-1) {
											neighbsurf[(nsurf%node)*6+4]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point)-1;
											k3 = sk + di;
											if (wall_on==0) {
												if(k3<0)   k3 += Nz;
												if(k3>=Nz) k3 -= Nz;
											}
											neighbsurf[(nsurf%node)*6+5]= 5*(i3 +(j3+k3*Ny)*Nx -myid*point);
										} else {
											k3 = sk + di;
											if (wall_on==0) {
												if(k3<0)   k3 += Nz;
												if(k3>=Nz) k3 -= Nz;
											}
											neighbsurf[(nsurf%node)*6+4]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point)-2;
											k3 = sk + 2*di;
											if (wall_on==0) {
												if(k3<0)   k3 += Nz;
												if(k3>=Nz) k3 -= Nz;
											}
											neighbsurf[(nsurf%node)*6+5]= 5*(i3 +(j3+k3*Ny)*Nx - myid*point);
										}
									}
								}
								nsurf++;
							}
						}
					}
				}
			}
		}
	}
	if (myid==0) printf("nsurf, nodes=%d %d\n",nsurf,nodes);
	if (Q_on!=0 && (patch_on!=0 || wall_on!=0) ) {
		MPI_Win_fence(0, winQsurf);
		MPI_Win_fence(0, winsurf);
		MPI_Win_fence(0, winneighbsurf);
	}
	if(flow_on!=0) MPI_Win_fence(0, winnf);
	MPI_Barrier(MPI_COMM_WORLD);	
}

