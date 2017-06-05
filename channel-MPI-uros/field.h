/*
 *  field.h
 *  
 *
 *  Created by Sirius on 12/12/16.
 *  Copyright 2015 home. All rights reserved.
 *
 * f_E: electric field
 * f_H: molecular field contributed by electric field
 * f_T: temperature field
 * f_ht: for evolving f_T 
 * f_Tsurf: temperature field at surface
 * f_htsurf: for evolving f_Tsurf
 */
#include "mpi.h"

extern int f_l, f_d;
extern double f_ea, f_Em, f_w0, f_lbd, f_betaE;
extern double f_k, f_zr;
extern double *f_E, *f_H, *f_T, *f_ht, *f_Tsurf, *f_htsurf;
extern MPI_Win winf_t, winf_e, winf_tsurf, winf_htsurf;

void f_allocate();
void f_deallocate();
void f_init();
void f_eq(double f_dt);
void f_evolveT(double *err, double f_dt);
void f_output3(int new);
