/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_coriolis.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "modify.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "fix_multisphere.h"  
#include "fix_relax_contacts.h"  
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;



/* ---------------------------------------------------------------------- */

FixCoriolis::FixCoriolis(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 13) error->all(FLERR,"Illegal fix gravity command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;


   if (strcmp(arg[3],"axis") == 0) {//on est sur lelement 4 et il regarde dans ce if si cest un vector si cest le cas il va ensuite verifier que la ligne aura bien 8 termes= on ajoute les 3 coordonnees du vecteur 
    xdir = force->numeric(FLERR,arg[4]);
    ydir = force->numeric(FLERR,arg[5]);
    zdir = force->numeric(FLERR,arg[6]);

    
  } else error->all(FLERR,"Illegal fix coriolis command");
  
  if (strcmp(arg[7],"point") == 0) {
    xpos = force->numeric(FLERR,arg[8]);
    ypos = force->numeric(FLERR,arg[9]);
    zpos = force->numeric(FLERR,arg[10]);

    
  } else error->all(FLERR,"Illegal fix coriolis command");
  
  if (strcmp(arg[11],"omega") == 0) {
    omega = force->numeric(FLERR,arg[12]);

    
  } else error->all(FLERR,"Illegal fix coriolis command");
  
  degree2rad = MY_PI/180.0;
  time_origin = update->ntimestep;

  eflag = 0;
  egrav = 0.0;

  fm = NULL; 
}

/* ---------------------------------------------------------------------- */

FixCoriolis::~FixCoriolis()
{

}

/* ---------------------------------------------------------------------- */

int FixCoriolis::setmask()
{
  int mask = 0;
  //mask |=POST_INTEGRATE;
  mask |= POST_FORCE;
  //mask |= END_OF_STEP;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCoriolis::init()
{

}


/* ---------------------------------------------------------------------- */

void FixCoriolis::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
  {
    post_force(vflag);
    //initial_integrate(vflag);
  }
}

/* ---------------------------------------------------------------------- */
/*
void FixCoriolis::initial_integrate(int vflag)
{
  // update gravity due to variables

  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double massone;


  eflag = 0;
  egrav = 0.0;


 if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (!fm || (fm && fm->belongs_to(i) < 0))) { 
        massone = rmass[i];
        f[i][0] += -2*massone*(ydir*v[i][2]-zdir*v[i][1])*omega;
        f[i][1] += -2*massone*(zdir*v[i][0]-xdir*v[i][2])*omega;
        f[i][2] += -2*massone*(xdir*v[i][1]-ydir*v[i][0])*omega;
        std::cout<<"endstep"<<f[i][0]<<std::endl;
	std::cout<<"endstep"<<f[i][1]<<std::endl;
	std::cout<<"endstep"<<f[i][2]<<std::endl;
      }
  } else {
	for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (!fm || (fm && fm->belongs_to(i) < 0))) { 
        massone = mass[type[i]];
        f[i][0] += -2*massone*(ydir*v[i][2]-zdir*v[i][1])*omega;
        f[i][1] += -2*massone*(zdir*v[i][0]-xdir*v[i][2])*omega;
        f[i][2] += -2*massone*(xdir*v[i][1]-ydir*v[i][0])*omega;
        //egrav -= massone * (xacc*x[i][0] + yacc*x[i][1] + zacc*x[i][2]);
      }
      }
      

}

*/
void FixCoriolis::post_force(int vflag)
{
  // update gravity due to variables

  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double massone;

  eflag = 0;
  egrav = 0.0;


 if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (!fm || (fm && fm->belongs_to(i) < 0))) { 
        massone = rmass[i];
        f[i][0] += -2*massone*(ydir*v[i][2]-zdir*v[i][1])*omega;
        f[i][1] += -2*massone*(zdir*v[i][0]-xdir*v[i][2])*omega;
        f[i][2] += -2*massone*(xdir*v[i][1]-ydir*v[i][0])*omega;
        /*std::cout<<"post_force"<<f[i][0]<<std::endl;
	std::cout<<"post_force"<<f[i][1]<<std::endl;
	std::cout<<"post_force"<<f[i][2]<<std::endl;*/
      }
  } else {
	for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (!fm || (fm && fm->belongs_to(i) < 0))) { 
        massone = mass[type[i]];
        f[i][0] += -2*massone*(ydir*v[i][2]-zdir*v[i][1])*omega;
        f[i][1] += -2*massone*(zdir*v[i][0]-xdir*v[i][2])*omega;
        f[i][2] += -2*massone*(xdir*v[i][1]-ydir*v[i][0])*omega;
        //egrav -= massone * (xacc*x[i][0] + yacc*x[i][1] + zacc*x[i][2]);
      }
      }
      

}
/*
void FixCoriolis::end_of_step()
{
  // update gravity due to variables

  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double massone;

  eflag = 0;
  egrav = 0.0;


 if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (!fm || (fm && fm->belongs_to(i) < 0))) { 
        massone = rmass[i];
        f[i][0] += -2*massone*(ydir*v[i][2]-zdir*v[i][1])*omega;
        f[i][1] += -2*massone*(zdir*v[i][0]-xdir*v[i][2])*omega;
        f[i][2] += -2*massone*(xdir*v[i][1]-ydir*v[i][0])*omega;
        /*std::cout<<"post_force"<<f[i][0]<<std::endl;
	std::cout<<"post_force"<<f[i][1]<<std::endl;
	std::cout<<"post_force"<<f[i][2]<<std::endl;*/
/*      }
  } else {
	for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (!fm || (fm && fm->belongs_to(i) < 0))) { 
        massone = mass[type[i]];
        f[i][0] += -2*massone*(ydir*v[i][2]-zdir*v[i][1])*omega;
        f[i][1] += -2*massone*(zdir*v[i][0]-xdir*v[i][2])*omega;
        f[i][2] += -2*massone*(xdir*v[i][1]-ydir*v[i][0])*omega;
        //egrav -= massone * (xacc*x[i][0] + yacc*x[i][1] + zacc*x[i][2]);
      }
      }
      

}
*/

/* ----------------------------------------------------------------------
   potential energy in gravity field
------------------------------------------------------------------------- */

double FixCoriolis::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(&egrav,&egrav_all,1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return egrav_all;
}

/* ---------------------------------------------------------------------- */


