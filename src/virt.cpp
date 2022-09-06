// ************************************************************************* //
// Virtual contributions to the cross section                                //
//                                                                           //
// By Benjamin Fuks - 18.07.2022                                             //
// ************************************************************************* //



// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
#include <complex>          // Complex numbers library                       //
// -------- Classes -------------------------------------------------------- //
#include "parameters.h"     // All the parameters of the calculation         //
#include "clooptools.h"     // Looptools                                     //
// ------------------------------------------------------------------------- //


// ************************************************************************* //
//  Virtual corrections.                                                     //
// ************************************************************************* //
double Virt_ll(const double s, const double t, const unsigned int flav1, const unsigned int flav2, Parameters *p)
{
   // PV integrals
   std::complex<double> c00[3], c1[3], c2[3], c12[3], c22[3];
   setmudim(p->muR2());
   for (int eps=0; eps<3; eps++)
     { setlambda(-eps); c00[eps]=C0i(cc00,0,0,s,0,0,0); c1[eps]=C0i(cc1,0,0,s,0,0,0); c2[eps]=C0i(cc2,0,0,s,0,0,0); c12[eps]=C0i(cc12,0,0,s,0,0,0); c22[eps]=C0i(cc22,0,0,s,0,0,0); }

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
     // eps^2/eps^2 pieces
     (LL+RR) * ( c00[2]*(12.*u2+4.*ut-6.*t2) + s2*((c12[2]+c2[2]+c22[2])*(2.*t-3.*u)+c1[2]*(t-u)) ) +
     (LR+RL) * ( c00[2]*(12.*t2+4.*ut-6.*u2) + s2*((c12[2]+c2[2]+c22[2])*(2.*u-3.*t)+c1[2]*(u-t)) );

   // result (including colour, flux and phase space)
   return 4.*p->alphas()/(3.*M_PI)*real(virt);
}

// ************************************************************************* //
//  Integrated dipole contributions.                                         //
// ************************************************************************* //
double Dipole_ll(const double s, const double t, const unsigned int flav1, const unsigned int flav2, Parameters *p)
{
   // Log and kinematics
   double u = -s-t, t2 = pow(t,2.), u2=pow(u,2.);
   std::complex<double> sz = s - p->cmz2();
   double lnmusq = std::log(p->muR2()/s);

   // couplings
   double LR =  norm(p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2)  / s + p->ZRll(p->flav1(),p->flav2())* p->ZLqq(flav1,flav2) / sz);
   double RL =  norm(p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2)  / s + p->ZLll(p->flav1(),p->flav2())* p->ZRqq(flav1,flav2) / sz);
   double LL =  norm(p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2)  / s + p->ZLll(p->flav1(),p->flav2())* p->ZLqq(flav1,flav2) / sz);
   double RR =  norm(p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2)  / s + p->ZRll(p->flav1(),p->flav2())* p->ZRqq(flav1,flav2) / sz);

   // Matrix element
   double dip =
     // eps^0 pieces
     (10. - pow(M_PI,2) + (3.+lnmusq)*lnmusq) * ((LR+RL)*t2 + (LL+RR)*u2) +
     // eps^-1 pieces
     (3. + 2.*lnmusq) * s * ((LR+RL)*(2.*t-u) + (LL+RR)*(2.*u-t)) +
     // eps^-2 pieces
     2. * (t2-u2) * (LR+RL-LL-RR);

   // result (including colour, flux and phase space)
   return 2.*p->alphas()/(3.*M_PI)*dip;
}

