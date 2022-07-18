// ************************************************************************* //
// NLO DY cross section                                                      //
//                                                                           //
// By Benjamin Fuks - 18.07.2022                                             //
// ************************************************************************* //



// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
#include <iostream>
// -------- Classes -------------------------------------------------------- //
#include "parameters.h"     // All the parameters of the calculation         //
#include "clooptools.h"     // Looptools                                     //
#include "dipoles.h"        // Dipole operators                              //
// -------- Functions ------------------------------------------------------ //
void Kinematics2to2(double &, double&, double&, double&, double&, double*x, ProcessConfig*);           //
void Kinematics2to2PK(double &, double&, double&, double*, double*, double*, double*x, ProcessConfig*);//
double Virt_ll(const double, const double, const unsigned int, const unsigned int, Parameters*);       //
double Dipole_ll(const double, const double, const unsigned int, const unsigned int, Parameters*);     //
double Born_Parton(const double &, const double &, const unsigned int &, Parameters*);                 //
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
        )* (Virt_ll(s,t,i,i,p) + Dipole_ll(s,t,i,i,p));
   }
   // disgma/dM
   if(p->M()!=-1) sig*=(2./p->M());

   // result (including colour, flux and phase space)
   return sig*jac/(48.*M_PI*pow(s,2.));
}



// ************************************************************************* //
//  Collinear reminder (P+K)                                                 //
// ************************************************************************* //
double PK(double *x, size_t dim, void *prm)
{
   // Conversion of type void* into Param*
   Parameters *p = (Parameters *)prm;

   // Kinematics
   double xa, xb, z, s[2], t[2], jac[3];
   Kinematics2to2PK(xa,xb,z,t,jac,s,x,p->Process());

   // Loop over the initial state quarks
   double sig=0.;
   for(unsigned int i=0; i<4; i++)
   {
      // Born partonic rates
      double born[2] = { Born_Parton(s[0],t[0],i,p), Born_Parton(s[1],t[1],i,p) };

      // P_qq and K_qq
      sig += p->alphas()/M_PI * (
          p->pdf()->xfxQ( i+1,xa,p->muF())*p->pdf()->xfxQ(-i-1,xb,p->muF())+
          p->pdf()->xfxQ(-i-1,xa,p->muF())*p->pdf()->xfxQ( i+1,xb,p->muF())
        )*(
         -Pqq_nodelta(z)                          * log(p->muF2()/s[1]) * jac[0] * born[1] -
          Pqq_delta(p->Process()->tauh()/(xa*xb)) * log(p->muF2()/s[0]) * jac[2] * born[0] +
          Pqq_plus(z)                             * log(p->muF2()/s[0]) * jac[1] * born[0] +
          Kqq_nodelta(z)                          * jac[0] * born[1] +
          Kqq_delta(p->Process()->tauh()/(xa*xb)) * jac[2] * born[0] -
          Kqq_plus(z)                             * jac[1] * born[0]
       );

      // P_gq and K_gq
      sig += p->alphas()/(2.*M_PI) * (
          p->pdf()->xfxQ(0,xa,p->muF()) * (p->pdf()->xfxQ(-i-1,xb,p->muF())+p->pdf()->xfxQ( i+1,xb,p->muF())) +
          p->pdf()->xfxQ(0,xb,p->muF()) * (p->pdf()->xfxQ(-i-1,xa,p->muF())+p->pdf()->xfxQ( i+1,xa,p->muF()))
        )*(
         -Pgq_nodelta(z)                          * log(p->muF2()/s[1]) * jac[0] * born[1] +
          Kgq_nodelta(z)                          * jac[0] * born[1]
       );
   }

   // Results
   return sig;
}




// ************************************************************************* //
//  Reals - dipole integrand.                                                //
// ************************************************************************* //
double Real(double *x, size_t dim, void *prm)
{
   // conversion of type void* into Param* and clearing looptools cache
   Parameters *p = (Parameters *)prm;

   // Kinematics
   double s=0., jac=0.;

  // Loop over the initial state quarks
   double sig=0.;

   // disgma/dM
   if(p->M()!=-1) sig*=(2./p->M());

   // result (including colour, flux and phase space)
   return 0.*sig*jac/(48.*M_PI*pow(s,2.));
}



