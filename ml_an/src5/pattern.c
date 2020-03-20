/*
 *  pattern.c
 *  
 *
 *  Created by Sirius on July/23/15.
 *  Copyright 2018. All rights reserved.
 *
 */

#include "pattern.h"
#include "main.h"

void patt_allocate()
{
    if (patt_on!=0) {
	    MPI_Win_allocate_shared(point*sizeof(real),1,MPI_INFO_NULL, shmcomm, &patt_c, &winpatt_c);
    }
}




void patt_deallocate()
{
    if (patt_on!=0) {
	    MPI_Win_fence(0, winpatt_c);
	    MPI_Win_free(&winpatt_c);
    }
}




void patt_init()
{
    FILE *param;
    int id, flag, top, i, j, k, norm, i0, j0, k0, ii, jj, ic;
    char junk[256], type;
    double rr, dx, dy, dz, sepx, sepy, sepz, cx, cy, cz, xx, yy;
    double px[1001], py[1001];

    if (patt_on!=0) {
//  initialize patt_c
        for (id=0; id<lpoint; id++) patt_c[id] = 0.;
	    MPI_Barrier(MPI_COMM_WORLD);	

//  read pattern.in
        param = fopen("pattern.in","r");

        if (param==NULL) {	// make sure param file exists
        	printf("input file pattern.in not exist\n");
        	exit(-1);
        } else if (myid==root) {
        	printf("Reading pattern.in ...\n");
        }
        
        flag = 1;		// skip info lines which start with #
        while ( flag==1 && fgets (junk , 256 , param) != NULL ) {
        	if(junk[0]=='#') {
        		flag=1;
        	} else {
        		flag=0;
        	}
        }

        do {
            sscanf(junk,"%c %d\n", &type, &top);
            if (type=='c') {                // circle
                double rc, width, pc;
                fgets (junk, 256, param);
                sscanf(junk,"%lf %lf %lf\n", &cx, &cy, &cz);
                fgets (junk , 256 , param);
                sscanf(junk,"r width %lf %lf\n", &rc, &width);
                for (k=0; k<Nz; k++) {
                    for (j=0; j<Ny; j++) {
                        for (i=0; i<Nx; i++) {
                            id = i + (j+k*Ny)*Nx;
                            if (id/point==myid) {
                                id = id%point;
                                dx = i - cx;
                                dy = j - cy;
                                dz = k - cz;
                                if (wall_x==0) {
                                    if (dx<-0.5*(double)Nx) dx+= Nx;
                                    else if (dx>0.5*(double)Nx) dx-=Nx;
                                }
                                if (wall_y==0) {
                                    if (dy<-0.5*(double)Ny) dy+= Ny;
                                    else if (dy>0.5*(double)Ny) dy-=Ny;
                                }
                                if (wall_z==0) {
                                    if (dz<-0.5*(double)Nz) dz+= Nz;
                                    else if (dz>0.5*(double)Nz) dz-=Nz;
                                }
                                if (Nx<=Ny && Nx<=Nz)       rr = sqrt(dy*dy+dz*dz);
                                else if (Ny<=Nx && Ny<=Nz)  rr = sqrt(dx*dx+dz*dz);
                                else                        rr = sqrt(dx*dx+dy*dy);
                                pc = 0.5 + 0.5*tanh((rr-rc)/width);
                                if (top<0) pc = 1.-pc;
                                if (pc>1e-7) patt_c[id] = pc;
                            }
                        }
                    }
                }
            } else if (type=='e') {                // ellipse
                double ra, rb, ira2, irb2, width, pc;
                fgets (junk, 256, param);
                sscanf(junk,"%lf %lf %lf\n", &cx, &cy, &cz);
                fgets (junk , 256 , param);
                sscanf(junk,"ra rb width %lf %lf %lf\n", &ra, &rb, &width);
                ira2 = 1./(ra*ra);
                irb2 = 1./(rb*rb);
                for (k=0; k<Nz; k++) {
                    for (j=0; j<Ny; j++) {
                        for (i=0; i<Nx; i++) {
                            id = i + (j+k*Ny)*Nx;
                            if (id/point==myid) {
                                id = id%point;
                                dx = i - cx;
                                dy = j - cy;
                                dz = k - cz;
                                if (wall_x==0) {
                                    if (dx<-0.5*(double)Nx) dx+= Nx;
                                    else if (dx>0.5*(double)Nx) dx-=Nx;
                                }
                                if (wall_y==0) {
                                    if (dy<-0.5*(double)Ny) dy+= Ny;
                                    else if (dy>0.5*(double)Ny) dy-=Ny;
                                }
                                if (wall_z==0) {
                                    if (dz<-0.5*(double)Nz) dz+= Nz;
                                    else if (dz>0.5*(double)Nz) dz-=Nz;
                                }
                                if (Nx<=Ny && Nx<=Nz)       rr = sqrt(dy*dy*ira2+dz*dz*irb2);
                                else if (Ny<=Nx && Ny<=Nz)  rr = sqrt(dx*dx*ira2+dz*dz*irb2);
                                else                        rr = sqrt(dx*dx*ira2+dy*dy*irb2);
                                pc = 0.5 + 0.5*tanh((rr-1.)/width);
                                if (top<0) pc = 1.-pc;
                                patt_c[id] = pc;
                            }
                        }
                    }
                }
            } else if (type=='d') {     // simple domain
                char dim;
                double cc, width, iwidth, dc, pc;
                fgets (junk, 256, param);
                sscanf(junk,"%c %lf %lf\n", &dim, &cc, &width);
                iwidth = 1./width;
                for (k=0; k<Nz; k++) {
                    for (j=0; j<Ny; j++) {
                        for (i=0; i<Nx; i++) {
                            id = i + (j+k*Ny)*Nx;
                            if (id/point==myid) {
                                id = id%point;
                                if (dim=='x')      dc = i - cc;
                                else if (dim=='y') dc = j - cc;
                                else               dc = k - cc;
                                pc = 0.5 + 0.5*tanh(dc*iwidth);
                                // close to lower boundary, the other boundary is at -0.5 or N-0.5
                                if (dim=='x')      dc = (double)i + 0.5;
                                else if (dim=='y') dc = (double)j + 0.5;
                                else               dc = (double)k + 0.5;
                                if (dc<0.5*cc) pc = 0.5 - 0.5*tanh(dc*iwidth);
                                // close to higher boundary
                                if (dim=='x')      dc = (double)Nx - 0.5 - i;
                                else if (dim=='y') dc = (double)Ny - 0.5 - j;
                                else               dc = (double)Nz - 0.5 - k;
                                if (dc<0.5*cc) pc = 0.5 + 0.5*tanh(dc*iwidth);
                                if (top<0) pc = 1. - pc;
                                patt_c[id] = pc;
                            }
                        }
                    }
                }
            } else if (type=='s') {     // square
                fgets (junk, 256, param);
                sscanf(junk,"%lf %lf %lf\n", &cx, &cy, &cz);
                fscanf(param, "%lf %lf %lf\n",&sepx,&sepy,&sepz);
                for (k=0; k<Nz; k++) {
                    for (j=0; j<Ny; j++) {
                        for (i=0; i<Nx; i++) {
                            id = i + (j+k*Ny)*Nx;
                            if (id/point==myid) {
                                id = id%point;
                                dx = i - cx;
                                dy = j - cy;
                                dz = k - cz;
                                if (wall_x==0) {
                                    if (dx<-0.5*(double)Nx) dx+= Nx;
                                    else if (dx>0.5*(double)Nx) dx-=Nx;
                                }
                                if (wall_y==0) {
                                    if (dy<-0.5*(double)Ny) dy+= Ny;
                                    else if (dy>0.5*(double)Ny) dy-=Ny;
                                }
                                if (wall_z==0) {
                                    if (dz<-0.5*(double)Nz) dz+= Nz;
                                    else if (dz>0.5*(double)Nz) dz-=Nz;
                                }
                                if (top>0 && fabs(dx*2.)<sepx && fabs(dy*2.)<sepy && fabs(dz*2.)<sepz) patt_c[id] = 1.;
                                else if (top<0 && (fabs(dx*2.)>=sepx || fabs(dy*2.)>=sepy || fabs(dz*2.)>=sepz)) patt_c[id] = 1.; 
                            }
                        }
                    }
                }
            } else if (type=='p') {     // polygon
                if (abs(top)>999) {
                    printf("polygon has too many corners\n");
                    exit(-1);
                }
                for (ii=0; ii<abs(top); ii++) fscanf(param, "%lf %lf\n", &px[ii], &py[ii]);
                norm = 2;               // normal of the quasi-2D system: 0: x, 1: y, 2: z.
                cx   = px[0];
                cy   = py[0];
                cz   = Nz/2;
                if (Nx<=Ny && Nx<=Nz) {
                    norm = 0;
                    cx   = Nx/2;
                    cy   = px[0];
                    cz   = py[0];
                }
                else if (Ny<=Nx && Ny<=Nz) {
                    norm = 1;
                    cx   = px[0];
                    cy   = Ny/2;
                    cz   = py[0];
                }
                for (k=0; k<Nz; k++) {
                    for (j=0; j<Ny; j++) {
                        for (i=0; i<Nx; i++) {
                            id = i + (j+k*Ny)*Nx;
                            if (id/point==myid) {
                                id = id%point;
//                               i0 = i;
//                               j0 = j;
//                               k0 = k;
//                               dx = i - cx;
//                               dy = j - cy;
//                               dz = k - cz;
//                               if (wall_x==0) {
//                                   if (dx<-0.5*(double)Nx) i0+= Nx;
//                                   else if (dx>0.5*(double)Nx) i0-=Nx;
//                               }
//                               if (wall_y==0) {
//                                   if (dy<-0.5*(double)Ny) j0+= Ny;
//                                   else if (dy>0.5*(double)Ny) j0-=Ny;
//                               }
//                               if (wall_z==0) {
//                                   if (dz<-0.5*(double)Nz) k0+= Nz;
//                                   else if (dz>0.5*(double)Nz) k0-=Nz;
//                               }

//                               xx = i0;
//                               yy = j0;
//                               if (norm!=2) {
//                                   yy = k0;
//                                   if (norm==0) xx = j0;
//                               }
                                xx = i;
                                yy = j;
                                if (norm!=2) {
                                    yy = k;
                                    if (norm==0) xx = j;
                                }
                                ic = 0;
                                for (ii=0, jj=abs(top)-1; ii<abs(top); jj=ii++) {
                                    if (((py[ii]<=yy && yy<py[jj]) || (py[jj]<=yy && yy<py[ii])) && xx<(px[jj]-px[ii])*(yy-py[ii])/(py[jj]-py[ii])+px[ii]) ic++; 
                                }

                                if (ic%2!=0 && top>0 || ic%2==0 && top<0) patt_c[id] = 1.;
                            }
                        }
                    }
                }
            } else if (type=='g') {     // smooth polygon
                double pc, width, dr1sq, dr2sq, dr3sq, drsqmin, lm;
                if (abs(top)>999) {
                    printf("polygon has too many corners\n");
                    exit(-1);
                }
                fscanf(param, "width %lf\n", &width);
                for (ii=0; ii<abs(top); ii++) fscanf(param, "%lf %lf\n", &px[ii], &py[ii]);
                px[ii] = px[0];
                py[ii] = py[0];
                if (ii!=abs(top)) printf("error in smooth polygon\n");
                norm = 2;               // normal of the quasi-2D system: 0: x, 1: y, 2: z.
                cx   = px[0];
                cy   = py[0];
                cz   = Nz/2;
                if (Nx<=Ny && Nx<=Nz) {
                    norm = 0;
                    cx   = Nx/2;
                    cy   = px[0];
                    cz   = py[0];
                }
                else if (Ny<=Nx && Ny<=Nz) {
                    norm = 1;
                    cx   = px[0];
                    cy   = Ny/2;
                    cz   = py[0];
                }
                for (k=0; k<Nz; k++) {
                    for (j=0; j<Ny; j++) {
                        for (i=0; i<Nx; i++) {
                            id = i + (j+k*Ny)*Nx;
                            if (id/point==myid) {
                                id = id%point;
                                xx = i;
                                yy = j;
                                if (norm!=2) {
                                    yy = k;
                                    if (norm==0) xx = j;
                                }
                                ic = 0;
                                for (ii=0, jj=abs(top)-1; ii<abs(top); jj=ii++) {
                                    if (((py[ii]<=yy && yy<py[jj]) || (py[jj]<=yy && yy<py[ii])) && xx<(px[jj]-px[ii])*(yy-py[ii])/(py[jj]-py[ii])+px[ii]) ic++; 
                                }

                                if (ic%2!=0 && top>0 || ic%2==0 && top<0) patt_c[id] = 1.;

                                dr2sq = (xx-px[0])*(xx-px[0]) + (yy-py[0])*(yy-py[0]);
                                drsqmin = dr2sq;
                                for (ii=0; ii<abs(top); ii++) {
                                    jj = ii+1;
                                    dr1sq = dr2sq;
                                    dr2sq = (xx-px[jj])*(xx-px[jj]) + (yy-py[jj])*(yy-py[jj]);

                                    lm = ((px[ii]-xx)*(px[ii]-px[jj])+(py[ii]-yy)*(py[ii]-py[jj]))/((py[ii]-py[jj])*(py[ii]-py[jj])+(px[ii]-px[jj])*(px[ii]-px[jj]));
                                    if (lm<0) {
                                        if (dr1sq<drsqmin) drsqmin=dr1sq;
                                    } else if (lm>1) {
                                        if (dr2sq<drsqmin) drsqmin=dr2sq;
                                    } else {
                                        dr3sq = (xx-(1.-lm)*px[ii]-lm*px[jj])*(xx-(1.-lm)*px[ii]-lm*px[jj])+(yy-(1.-lm)*py[ii]-lm*py[jj])*(yy-(1.-lm)*py[ii]-lm*py[jj]);
                                        if (dr3sq<drsqmin) drsqmin=dr3sq;
                                    }
                                }
                                if (patt_c[id]<0.5) pc = 0.5+0.5*tanh(-sqrt(drsqmin)/width);
                                else                pc = 0.5+0.5*tanh( sqrt(drsqmin)/width);
                                patt_c[id] = pc;
                            }
                        }
                    }
                }
            } else if (type=='f') {     // flip passive & active region
                for (id=0; id<lpoint; id++) patt_c[id] = 1.-patt_c[id];
            }
	        MPI_Win_fence(0,winpatt_c);
	        MPI_Barrier(MPI_COMM_WORLD);	
        } while ( fgets (junk , 256 , param) != NULL ); 

        fclose(param);
	    MPI_Barrier(MPI_COMM_WORLD);	
    }
}
