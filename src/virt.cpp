// ************************************************************************* //
// Virtual contributions to the cross section                                //
//                                                                           //
// By Benjamin Fuks - 08.07.2022                                             //
// ************************************************************************* //



// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
#include <complex>          // Complex numbers library                       //
#include <iostream>
// -------- Classes -------------------------------------------------------- //
#include "parameters.h"     // All the parameters of the calculation         //
#include "clooptools.h"     // Looptools                                     //
// ------------------------------------------------------------------------- //


// ************************************************************************* //
//  Virtual corrections.                                                     //
// ************************************************************************* //
double virt_ll(const double s, const double t, unsigned int flav1, unsigned int flav2, Parameters *p)
{
   // PV integrals
   std::complex<double> c00[3], c1[3], c2[3], c12[3], c22[3];

   setmudim(p->mz2());
   setlambda( 0); c00[0] = C0i(cc00, 0., 0., s, 0., 0., 0.); c1[0]  = C0i(cc1,  0., 0., s, 0., 0., 0.); c2[0]  = C0i(cc2,  0., 0., s, 0., 0., 0.); c12[0] = C0i(cc12, 0., 0., s, 0., 0., 0.); c22[0] = C0i(cc22, 0., 0., s, 0., 0., 0.);
   setlambda(-1); c00[1] = C0i(cc00, 0., 0., s, 0., 0., 0.); c1[1]  = C0i(cc1,  0., 0., s, 0., 0., 0.); c2[1]  = C0i(cc2,  0., 0., s, 0., 0., 0.); c12[1] = C0i(cc12, 0., 0., s, 0., 0., 0.); c22[1] = C0i(cc22, 0., 0., s, 0., 0., 0.);
   setlambda(-2); c00[2] = C0i(cc00, 0., 0., s, 0., 0., 0.); c1[2]  = C0i(cc1,  0., 0., s, 0., 0., 0.); c2[2]  = C0i(cc2,  0., 0., s, 0., 0., 0.); c12[2] = C0i(cc12, 0., 0., s, 0., 0., 0.); c22[2] = C0i(cc22, 0., 0., s, 0., 0., 0.);

   // mandelstam and propagators
   double u = -s-t, t2 = pow(t,2.), u2=pow(u,2.), s2=pow(s,2.), ut=u*t;
   std::complex<double> sz = s - p->cmz2();

   // Couplings
   double LR = norm( p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2) / s + p->ZRll(p->flav1(),p->flav2())* p->ZLqq(flav1,flav2) / sz );
   double RL = norm( p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2) / s + p->ZLll(p->flav1(),p->flav2())* p->ZRqq(flav1,flav2) / sz );
   double LL = norm( p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2) / s + p->ZLll(p->flav1(),p->flav2())* p->ZLqq(flav1,flav2) / sz );
   double RR = norm( p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2) / s + p->ZRll(p->flav1(),p->flav2())* p->ZRqq(flav1,flav2) / sz );

   // Matrix element
   std::complex<double> virt =
     // eps^0 pieces
     2. * ( (LL+RR)*u2 + (LR+RL)*t2) * (c00[0] + (c1[0]+c12[0]+c2[0]+c22[0])*s/2.) +
     // eps^1/eps pieces
     (LL+RR) * ( 2.*c00[1]*(t2-ut-4.*u2) + s*((c12[1]+c2[1]+c22[1])*(t2-ut-3.*u2)+c1[1]*(t2-ut-2.*u2)) ) +
     (LR+RL) * ( 2.*c00[1]*(u2-ut-4.*t2) + s*((c12[1]+c2[1]+c22[1])*(u2-ut-3.*t2)+c1[1]*(u2-ut-2.*t2)) ) +
     // eps^1/eps pieces
     (LL+RR) * ( c00[2]*(12.*u2+4.*ut-6.*t2) + s2*((c12[2]+c2[2]+c22[2])*(2.*t-3.*u)+c1[2]*(t-u)) ) +
     (LR+RL) * ( c00[2]*(12.*t2+4.*ut-6.*u2) + s2*((c12[2]+c2[2]+c22[2])*(2.*u-3.*t)+c1[2]*(u-t)) );

   // result (including colour, flux and phase space)
   return 4.*p->alphas()/(3.*M_PI)*real(virt);
}

// ************************************************************************* //
//  Integrated dipole contributions.                                         //
// ************************************************************************* //
double dipole_ll(const double s, const double t, Parameters *p)
{
   // result (including colour, flux and phase space)
   return 0.;
}

