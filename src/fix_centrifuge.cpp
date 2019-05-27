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
#include "fix_centrifuge.h"
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

FixCentrifuge::FixCentrifuge(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 13) error->all(FLERR,"Illegal fix centrifugal command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;


   if (strcmp(arg[3],"axis") == 0) {//on est sur lelement 4 et il regarde dans ce if si cest un vector si cest le cas il va ensuite verifier que la ligne aura bien 8 termes= on ajoute les 3 coordonnees du vecteur 
    xdir = force->numeric(FLERR,arg[4]);
    ydir = force->numeric(FLERR,arg[5]);
    zdir = force->numeric(FLERR,arg[6]);

    
  } else error->all(FLERR,"Illegal fix centrifugal command");
  
  if (strcmp(arg[7],"point") == 0) {
    xpos = force->numeric(FLERR,arg[8]);
    ypos = force->numeric(FLERR,arg[9]);
    zpos = force->numeric(FLERR,arg[10]);

    
  } else error->all(FLERR,"Illegal fix centrifugal command");
  
  if (strcmp(arg[11],"omega") == 0) {
    omega = force->numeric(FLERR,arg[12]);

    
  } else error->all(FLERR,"Illegal fix centrifugal command");
  
  degree2rad = MY_PI/180.0;
  time_origin = update->ntimestep;

  eflag = 0;
  

  fm = NULL; 
}

/* ---------------------------------------------------------------------- */

FixCentrifuge::~FixCentrifuge()
{

}

/* ---------------------------------------------------------------------- */

int FixCentrifuge::setmask()
{
  int mask = 0;
  //mask |= END_OF_STEP;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCentrifuge::init()
{

}


/* ---------------------------------------------------------------------- */

void FixCentrifuge::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")){
    post_force(vflag);
    //initial_integrate(vflag);
    }
}

/* ---------------------------------------------------------------------- */


void FixCentrifuge::post_force(int vflag)
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
  double t;
  double xH[3];
  double direction[3];
  double distance;
  
  //Determination des coordonnee du point de projection orthogonal (on passe par l'equation parametrique de l'axe de rotation)
 
  
  

  eflag = 0;


 if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (!fm || (fm && fm->belongs_to(i) < 0))) { 
        massone = rmass[i];
        t=-((xpos-x[i][0])*xdir+(ypos-x[i][1])*ydir+(zpos-x[i][2])*zdir)/(xdir*xdir+ydir*ydir+zdir*zdir);
  	xH[0]=xdir*t+xpos;
  	xH[1]=ydir*t+ypos;
  	xH[2]=zdir*t+zpos;
  	direction[0]=x[i][0]-xH[0];
  	direction[1]=x[i][1]-xH[1];
  	direction[2]=x[i][2]-xH[2];
  	distance=sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
        f[i][0] += massone*omega*omega*direction[0];
        f[i][1] += massone*omega*omega*direction[1];
        f[i][2] += massone*omega*omega*direction[2];
        //std::cout<<"distance"<<distance<<std::endl;
      }
  } else {
	for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (!fm || (fm && fm->belongs_to(i) < 0))) { 
        massone = mass[type[i]];
        t=-((xpos-x[i][0])*xdir+(ypos-x[i][1])*ydir+(zpos-x[i][2])*zdir)/(xdir*xdir+ydir*ydir+zdir*zdir);
  	xH[0]=xdir*t+xpos;
  	xH[1]=ydir*t+ypos;
  	xH[2]=zdir*t+zpos;
  	direction[0]=x[i][0]-xH[0];
  	direction[1]=x[i][1]-xH[1];
  	direction[2]=x[i][2]-xH[2];
  	distance=sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
        f[i][0] += massone*omega*omega*direction[0];
        f[i][1] += massone*omega*omega*direction[1];
        f[i][2] += massone*omega*omega*direction[2];
      }
      }
      
      

}





/* ----------------------------------------------------------------------
   potential energy in gravity field
------------------------------------------------------------------------- */

double FixCentrifuge::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(&egrav,&egrav_all,1,MPI_DOUBLE,MPI_SUM,world);
    eflag = 1;
  }
  return egrav_all;
}

/* ---------------------------------------------------------------------- */


