/*
 *  particle.h
 *  
 *
 *  Created by Sirius on 3/3/15.
 *  Copyright 2015 home. All rights reserved.
 *
 */

extern double Lx, Ly, Lz, hLx, hLy, hLz;
extern double *p_m, *p_I, *p_rad, *p_pos, *p_vel, *p_acc, *p_qua, *p_angv, *p_torq;
extern double *p_Wtype, *p_W, *ubounce;
extern int MAXPpoint;
extern int *plink, *link_plist;

void p_allocate();
void p_deallocate();
void p_init();
void p_iden();
void p_update();
int onsurface(double x, double y, double z, int ipar);