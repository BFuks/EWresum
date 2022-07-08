// ************************************************************************* //
// NLO DY cross section                                                      //
//                                                                           //
// By Benjamin Fuks - 08.07.2022                                             //
// ************************************************************************* //



// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Classes -------------------------------------------------------- //
#include "parameters.h"     // All the parameters of the calculation         //
#include "clooptools.h"     // Looptools                                     //
// -------- Functions ------------------------------------------------------ //
void Kinematics2to2(double &, double&, double&, double&, double&, double*x, ProcessConfig*); //
double virt_ll(const double, const double, unsigned int, unsigned int, Parameters*);         //
double dipole_ll(const double, const double, Parameters*);                   //
// ------------------------------------------------------------------------- //



// ************************************************************************* //
//  Virtual + dipole integrand.                                              //
// ************************************************************************* //
double Virt(double *x, size_t dim, void *prm)
{
   // conversion of type void* into Param* and clearing looptools cache
   Parameters *p = (Parameters *)prm;
   clearcache();

   // Kinematics 
   double xa, xb, s, t, jac;
   Kinematics2to2(xa,xb,t,jac,s,x,p->Process());

   // Loop over the initial state quarks
   double sig=0.;
   for(unsigned int i=0; i<4; i++)
   {
      sig+=(
          p->pdf()->xfxQ( i+1,xa,p->muF())*p->pdf()->xfxQ(-i-1,xb,p->muF())+
          p->pdf()->xfxQ(-i-1,xa,p->muF())*p->pdf()->xfxQ( i+1,xb,p->muF())
        )* (virt_ll(s,t,i,i,p) + dipole_ll(s,t,p));
   }

   // disgma/dM
   if(p->M()!=-1) sig*=(2./p->M());

   // result (including colour, flux and phase space)
   return sig*jac/(48.*M_PI*pow(s,2.));
}

