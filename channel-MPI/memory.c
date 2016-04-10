/*
 *  memory.c
 *  
 *
 *  Created by Sirius on 2/21/15.
 *  Copyright 2015 home. All rights reserved.
 *
 */
#include "main.h"

void allocate()
{
        if (myid==root) printf("Allocating..\n");

	MPI_Win_allocate_shared(5*point*sizeof(real),1,MPI_INFO_NULL, shmcomm, &Q,   &winq);
	MPI_Win_allocate_shared(  point*sizeof(real),1,MPI_INFO_NULL, shmcomm, &Rho, &winr);
	MPI_Win_allocate_shared(3*point*sizeof(real),1,MPI_INFO_NULL, shmcomm, &u,   &winu);

	if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) {
		MPI_Win_allocate_shared( 5*node*sizeof(real),1,MPI_INFO_NULL, shmcomm, &Qsurf, &winQsurf);
		MPI_Win_allocate_shared( 5*node*sizeof(real),1,MPI_INFO_NULL, shmcomm, &Hsurf, &winHsurf);
		MPI_Win_allocate_shared( 6*node*sizeof(int), 1,MPI_INFO_NULL, shmcomm, &neighbsurf, &winneighbsurf);
		MPI_Win_allocate_shared(10*node*sizeof(real),1,MPI_INFO_NULL, shmcomm, &surf, &winsurf);
	}

	if (flow_on!=0 && Q_on!=0) {
		MPI_Win_allocate_shared( 9*point*sizeof(real),1,MPI_INFO_NULL, shmcomm, &sigma_q, &wins);
	}
	if (flow_on!=0) {
		MPI_Win_allocate_shared(15*point*sizeof(real),1,MPI_INFO_NULL, shmcomm, &f, &winf);
                MPI_Win_allocate_shared(15*point*sizeof(real),1,MPI_INFO_NULL, shmcomm, &f2,&winf2);
		MPI_Win_allocate_shared(15*point*sizeof(real),1,MPI_INFO_NULL, shmcomm, &p, &winp);
		MPI_Win_allocate_shared(15*point*sizeof(int), 1,MPI_INFO_NULL, shmcomm, &nextf, &winnf);
	}

	info= malloc(point*sizeof(int));
	
	if (Q_on!=0) {
		neighb  = malloc( 6*point*sizeof(int));
		H       = malloc( 5*point*sizeof(real));
	}
	
	if (flow_on!=0) {
		Cf      = malloc(15*point*sizeof(real));
	}

	if (flow_on!=0 && Q_on!=0) {
		sigma_p = malloc( 3*point*sizeof(real));
		W       = malloc( 9*point*sizeof(real));
	}

	MPI_Barrier(MPI_COMM_WORLD);
}


void deallocate()
{
	if (myid==root) printf("Deallocating..\n");

        MPI_Win_fence(0, winq);
        MPI_Win_fence(0, winr);
        MPI_Win_fence(0, winu);
        MPI_Win_free(&winq);
        MPI_Win_free(&winr);
        MPI_Win_free(&winu);

	if (wall_x!=0 || wall_y!=0 || wall_z!=0 || npar>0) {
		MPI_Win_fence(0, winQsurf);
		MPI_Win_fence(0, winHsurf);
		MPI_Win_fence(0, winneighbsurf);
		MPI_Win_fence(0, winsurf);
                MPI_Win_free(&winQsurf);
                MPI_Win_free(&winHsurf);
                MPI_Win_free(&winneighbsurf);
                MPI_Win_free(&winsurf);
        }
	if(Q_on!=0 && flow_on!=0) {
		MPI_Win_fence(0, wins);
		MPI_Win_free(&wins);
	}

        if(flow_on!=0){
		MPI_Win_fence(0, winf);
		MPI_Win_fence(0, winf2);
		MPI_Win_fence(0, winp);
		MPI_Win_fence(0, winnf);
                MPI_Win_free(&winf);
                MPI_Win_free(&winf2);
                MPI_Win_free(&winp);
		MPI_Win_free(&winnf);
        }

	free(info);

	if (Q_on!=0) {		
		free(H);
		free(neighb);
	}
	
	if (flow_on!=0) {
		free(Cf);
	}

	if (flow_on!=0 && Q_on!=0) {
		free(sigma_p);
		free(W);
	}

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Comm_free(&shmcomm);
        MPI_Finalize();
}
