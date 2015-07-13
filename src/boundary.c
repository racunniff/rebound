/**
 * @file 	boundary.c
 * @brief 	Implementation of all boundary conditions. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The code supports different boundary conditions.
 * 
 * 
 * @section LICENSE
 * Copyright (c) 2015 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "integrator.h"
#include "rebound.h"
#include "boundary.h"
#include "tree.h"

void boundary_check(struct Rebound* const r){
	struct reb_particle* const particles = r->particles;
	const int N = r->N;
	const double boxsize_x = r->boxsize_x;
	const double boxsize_y = r->boxsize_y;
	const double boxsize_z = r->boxsize_z;
	switch(r->boundary){
		case RB_BT_OPEN:
			for (int i=0;i<N;i++){ // run through loop backwards so we don't have to recheck same index after removing
				int removep = 0;
				if(particles[i].x>boxsize_x/2.){
					removep = 1;
				}
				if(particles[i].x<-boxsize_x/2.){
					removep = 1;
				}
				if(particles[i].y>boxsize_y/2.){
					removep = 1;
				}
				if(particles[i].y<-boxsize_y/2.){
					removep = 1;
				}
				if(particles[i].z>boxsize_z/2.){
					removep = 1;
				}
				if(particles[i].z<-boxsize_z/2.){
					removep = 1;
				}
				if (removep==1){
					if (r->tree_root==NULL){
						particles_remove(r, i,0); // keepSorted=0 by default in C version
						i--; // need to recheck the particle that replaced the removed one
					}else{
						fprintf(stderr,"\n\033[1mError!\033[0m Cannot remove particle from tree.");
						exit(EXIT_FAILURE);
					}
				}
			}
			break;
		case RB_BT_SHEAR:
		{
			// The offset of ghostcell is time dependent.
			const double OMEGA = r->ri_sei.OMEGA;
			const double offsetp1 = -fmod(-1.5*OMEGA*boxsize_x*r->t+boxsize_y/2.,boxsize_y)-boxsize_y/2.; 
			const double offsetm1 = -fmod( 1.5*OMEGA*boxsize_x*r->t-boxsize_y/2.,boxsize_y)+boxsize_y/2.; 
			struct reb_particle* const particles = r->particles;
#pragma omp parallel for schedule(guided)
			for (int i=0;i<N;i++){
				// Radial
				while(particles[i].x>boxsize_x/2.){
					particles[i].x -= boxsize_x;
					particles[i].y += offsetp1;
					particles[i].vy += 3./2.*OMEGA*boxsize_x;
				}
				while(particles[i].x<-boxsize_x/2.){
					particles[i].x += boxsize_x;
					particles[i].y += offsetm1;
					particles[i].vy -= 3./2.*OMEGA*boxsize_x;
				}
				// Azimuthal
				while(particles[i].y>boxsize_y/2.){
					particles[i].y -= boxsize_y;
				}
				while(particles[i].y<-boxsize_y/2.){
					particles[i].y += boxsize_y;
				}
				// Vertical (there should be no boundary, but periodic makes life easier)
				while(particles[i].z>boxsize_z/2.){
					particles[i].z -= boxsize_z;
				}
				while(particles[i].z<-boxsize_z/2.){
					particles[i].z += boxsize_z;
				}
			}
		}
		break;
		case RB_BT_PERIODIC:
#pragma omp parallel for schedule(guided)
			for (int i=0;i<N;i++){
				while(particles[i].x>boxsize_x/2.){
					particles[i].x -= boxsize_x;
				}
				while(particles[i].x<-boxsize_x/2.){
					particles[i].x += boxsize_x;
				}
				while(particles[i].y>boxsize_y/2.){
					particles[i].y -= boxsize_y;
				}
				while(particles[i].y<-boxsize_y/2.){
					particles[i].y += boxsize_y;
				}
				while(particles[i].z>boxsize_z/2.){
					particles[i].z -= boxsize_z;
				}
				while(particles[i].z<-boxsize_z/2.){
					particles[i].z += boxsize_z;
				}
			}
		break;
		default:
		break;
	}
}

const static struct Ghostbox nan_ghostbox = {.shiftx = 0, .shifty = 0, .shiftz = 0, .shiftvx = 0, .shiftvy = 0, .shiftvz = 0};

struct Ghostbox boundary_get_ghostbox(struct Rebound* const r, int i, int j, int k){
	switch(r->boundary){
		case RB_BT_OPEN:
		{
			struct Ghostbox gb;
			gb.shiftx = r->boxsize_x*(double)i;
			gb.shifty = r->boxsize_y*(double)j;
			gb.shiftz = r->boxsize_z*(double)k;
			gb.shiftvx = 0;
			gb.shiftvy = 0;
			gb.shiftvz = 0;
			return gb;
		}
		case RB_BT_SHEAR:
		{
			const double OMEGA = r->ri_sei.OMEGA;
			struct Ghostbox gb;
			// Ghostboxes habe a finite velocity.
			gb.shiftvx = 0.;
			gb.shiftvy = -1.5*(double)i*OMEGA*r->boxsize_x;
			gb.shiftvz = 0.;
			// The shift in the y direction is time dependent. 
			double shift;
			if (i==0){
				shift = -fmod(gb.shiftvy*r->t,r->boxsize_y); 
			}else{
				if (i>0){
					shift = -fmod(gb.shiftvy*r->t-r->boxsize_y/2.,r->boxsize_y)-r->boxsize_y/2.; 
				}else{
					shift = -fmod(gb.shiftvy*r->t+r->boxsize_y/2.,r->boxsize_y)+r->boxsize_y/2.; 
				}	
			}
			gb.shiftx = r->boxsize_x*(double)i;
			gb.shifty = r->boxsize_y*(double)j-shift;
			gb.shiftz = r->boxsize_z*(double)k;
			return gb;
		}
		case RB_BT_PERIODIC:
		{
			struct Ghostbox gb;
			gb.shiftx = r->boxsize_x*(double)i;
			gb.shifty = r->boxsize_y*(double)j;
			gb.shiftz = r->boxsize_z*(double)k;
			gb.shiftvx = 0;
			gb.shiftvy = 0;
			gb.shiftvz = 0;
			return gb;
		}
		default:
			return nan_ghostbox;
	}
}

/**
 * Checks if a given particle is within the computational domain.
 * @param p reb_particle to be checked.
 * @return Return value is 1 if particle is inside the box and 0 otherwise.
 */
int boundary_particle_is_in_box(const struct Rebound* const r, struct reb_particle p){
	switch(r->boundary){
		case RB_BT_OPEN:
		case RB_BT_SHEAR:
			if(p.x>r->boxsize_x/2.){
				return 0;
			}
			if(p.x<-r->boxsize_x/2.){
				return 0;
			}
			if(p.y>r->boxsize_y/2.){
				return 0;
			}
			if(p.y<-r->boxsize_y/2.){
				return 0;
			}
			if(p.z>r->boxsize_z/2.){
				return 0;
			}
			if(p.z<-r->boxsize_z/2.){
				return 0;
			}
			return 1;
		default:
			return 1;
	}
}

