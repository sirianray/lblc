/*
This source file is to list all the lattice points that are within a helix


*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
	int a[41][41][41];
	int i, j, k, n;
	double nn = 2.65, b = 5.7, r = 4.01, th = 0.99, ph0=2.0;
	double xc = 20.5, yc = 20.5, zc = 20.5;
	double t, twopi = 6.2831853071795862;
	double x, y, z, dx, dy, dz, dr2;
	
	for (i=0; i<41; i++) {
		for (j=0; j<41; j++) {
			for (k=0; k<41; k++) {
				a[i][j][k]=0;
			}
		}
	}

	for (n=-2500; n<2500; n++) {
		t = 0.001 * n;
		x = xc + r*cos(twopi * t + ph0 );
		z = zc - r*sin(twopi * t + ph0 );
		y = yc + b * t;
		for (i=(int)(x-th-1.0); i<=(int)(x+th+1.0); i++) {
			for (j=(int)(y-th-1.0); j<=(int)(y+th+1.0); j++) {
				for (k=(int)(z-th-1.0); k<=(int)(z+th+1.0); k++) {
					dx = x - i;
					dy = y - j;
					dz = z - k;
					dr2= dx*dx + dy*dy +dz*dz;
					if(dr2<th*th){
						a[i][j][k]=1;
//						if(i==16 && j==20 & k==11)printf("%f %f %f %f 99768\n",x,y,z,t);
					}
				}
			}
		}
	}

        for (i=0; i<41; i++) {
                for (j=0; j<41; j++) {
                        for (k=0; k<41; k++) {
                                if(a[i][j][k]==1)printf("%d %d %d\n",i,j,k);
                        }
                }
        }



	return 0;
}
