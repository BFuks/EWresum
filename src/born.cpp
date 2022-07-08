// ************************************************************************* //
// Born cross section                                                        //
//                                                                           //
// By Benjamin Fuks - 08.07.2022                                             //
// ************************************************************************* //



// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers ------------------------------------- //
#include <complex>          // Complex numbers library                       //
// -------- Classes -------------------------------------------------------- //
#include "parameters.h"     // All the parameters of the calculation         //
// -------- Functions ------------------------------------------------------ //
void Kinematics2to2(double &, double&, double&, double&, double&, double*x, ProcessConfig*);
// ------------------------------------------------------------------------- //



// ************************************************************************* //
//  Born integrand .                                                         //
// ************************************************************************* //
double Born(double *x, size_t dim, void *prm)
{
   // conversion of type void* into Param*
   Parameters *p = (Parameters *)prm;

   // Kinematics 
   double xa, xb, s, t, jac;
   Kinematics2to2(xa,xb,t,jac,s,x,p->Process());
   double u  = -s-t;
   std::complex<double> sz = s - p->cmz2();

   // Loop over the initial state quarks
   double sig=0.;
   for(unsigned int i=0; i<4; i++)
   {
      double LR =  norm(p->All(p->flav1(),p->flav2()) * p->Aqq(i,i)  / s +
                      p->ZRll(p->flav1(),p->flav2())* p->ZLqq(i,i) / sz);
      double RL =  norm(p->All(p->flav1(),p->flav2()) * p->Aqq(i,i)  / s +
                      p->ZLll(p->flav1(),p->flav2())* p->ZRqq(i,i) / sz);
      double LL =  norm(p->All(p->flav1(),p->flav2()) * p->Aqq(i,i)  / s +
                      p->ZLll(p->flav1(),p->flav2())* p->ZLqq(i,i) / sz);
      double RR =  norm(p->All(p->flav1(),p->flav2()) * p->Aqq(i,i)  / s +
                      p->ZRll(p->flav1(),p->flav2())* p->ZRqq(i,i) / sz);
      sig+=(
          p->pdf()->xfxQ( i+1,xa,p->muF())*p->pdf()->xfxQ(-i-1,xb,p->muF())+
          p->pdf()->xfxQ(-i-1,xa,p->muF())*p->pdf()->xfxQ( i+1,xb,p->muF())
        )* (pow(t,2.) * (LR + RL) + pow(u,2.) * (LL + RR));
   }

   // disgma/dM
   if(p->M()!=-1) sig*=(2./p->M());

   // result (including colour, flux and phase space)
   return sig/(48.*M_PI*pow(s,2.))*jac;
}

