// ************************************************************************* //
// Real emission contributions to the cross section                          //
//                                                                           //
// By Benjamin Fuks - 06.09.2022                                             //
// ************************************************************************* //



// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
#include <complex>          // Complex numbers library                       //
// -------- Classes -------------------------------------------------------- //
#include "parameters.h"     // All the parameters of the calculation         //
// -------- Functions ------------------------------------------------------ //
double Born_Parton(const double&, const double &, const unsigned int &, Parameters *);//
// ------------------------------------------------------------------------- //


// ************************************************************************* //
//  Real gluon emission contributions                                        //
// ************************************************************************* //
double Real_Gluon_ll(const double s, const double M2, const double pt2, const double cth, const double phi, const unsigned int flav1, const unsigned int flav2, Parameters *p)
{
   // Kinematics
   double pap3[2] = { (s - M2 + sqrt(pow(s-M2,2) - 4.*s*pt2)) / 4., (s - M2 - sqrt(pow(s-M2,2) - 4.*s*pt2)) / 4.};
   double pap1[2], pbp1[2];
   for (unsigned int i=0;i<2;i++)
   {
     pap1[i] = (s - 2.*pap3[i]) * (1.-cth) / 4.;
     double cA = 1. - 2. * s * M2 / ( (M2 + 2.*pap3[i]) * (s - 2.*pap3[i]));
     pbp1[i] = .25 * (M2+ 2.*pap3[i]) * (1.-cA*cth-sqrt(1.-cA*cA)*sqrt(1.-cth*cth)*cos(phi));
   }
   std::complex<double> M2z = M2 - p->cmz2();

   // Couplings
   double LR = norm( p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2) / M2 + p->ZRll(p->flav1(),p->flav2())* p->ZLqq(flav1,flav2) / M2z );
   double RL = norm( p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2) / M2 + p->ZLll(p->flav1(),p->flav2())* p->ZRqq(flav1,flav2) / M2z );
   double LL = norm( p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2) / M2 + p->ZLll(p->flav1(),p->flav2())* p->ZLqq(flav1,flav2) / M2z );
   double RR = norm( p->All(p->flav1(),p->flav2()) * p->Aqq(flav1,flav2) / M2 + p->ZRll(p->flav1(),p->flav2())* p->ZRqq(flav1,flav2) / M2z );

   // matrix element, summing over two rapidities
   double realgg = 0.;
   for (unsigned int i=0; i<2; i++)
     realgg += -1. / pap3[i] / (M2-s+2.*pap3[i]) * (
       (LL+RR) * ( pow(s -2.*pap1[i]-2.*pap3[i],2)+4.*pow(pbp1[i],2)) +
       (LR+RL) * ( pow(M2-2.*pbp1[i]+2.*pap3[i],2)+4.*pow(pap1[i],2)) );

   // Output
   return 2.*p->alphas()/M_PI*M2*realgg;
}


double Dipole_Gluon_ll(const double s, const double M2, const double pt2, const double cth, const double phi, const unsigned int flav1, const unsigned int flav2, Parameters *p)
{
   // General dot products
   double pap3[2] = { (s - M2 + sqrt(pow(s-M2,2) - 4.*s*pt2)) / 4., (s - M2 - sqrt(pow(s-M2,2) - 4.*s*pt2)) / 4.};
   double pap1[2], pbp1[2], p1p3[2];

   // Leg 1 - other dot products + rescaling of pap1 (dipole kinematics)
   for (unsigned int i=0;i<2;i++)
   {
     pap1[i] = (s - 2.*pap3[i]) * (1.-cth) / 4.;
     const double cA = 1. - 2. * s * M2 / ( (M2 + 2.*pap3[i]) * (s - 2.*pap3[i]));
     pbp1[i] = .25 * (M2+ 2.*pap3[i]) * (1.-cA*cth-sqrt(1.-cA*cA)*sqrt(1.-cth*cth)*cos(phi));
     p1p3[i] = pap1[i]+pbp1[i]-M2/2.;

     // Rescaling of pap1 to the dipole kinematics
     pap1[i] = (2.*M2*(M2*pap1[i] + pap3[i]*(pap1[i] + pbp1[i])) + M2*(M2 + pap3[i] - 2.*pbp1[i])*s)/(2.*(2.*M2 + (1.0 - M2/s) * pap3[i])*s);
   }

   // Leg 1 - Matrix element, summing over two rapidities
   double dipg = 0.;
   for (unsigned int i=0; i<2; i++) dipg += Born_Parton(M2,-2.*pap1[i],flav1,p)/pap3[i];

   // Leg 2 - other dot products + rescaling of pap1 (dipole kinematics)
   for (unsigned int i=0;i<2;i++)
   {
     // reset of pap1
     pap1[i] = (s - 2.*pap3[i]) * (1.-cth) / 4.;

     // Rescaling pbp1 to the dipole kinematics
     const double den = 2.*M2 + (1.-M2/s)*((s-M2)/2.-pap3[i]);
     pbp1[i] = M2*(-4.*pap3[i]*(pap1[i] + pbp1[i]) - 2*(pap1[i] + pap3[i] - pbp1[i])*s + pow(s,2) + M2*(-2.*pap1[i] + 2*pbp1[i] + s))/(4.*den*s);

     // Rescaling pap1 to the dipole kinematics
     pap1[i] = (M2*(-(M2*p1p3[i]) + 2*p1p3[i]*pap3[i] + 2*M2*(pap1[i] + pap3[i]) - p1p3[i]*s + 2*pap1[i]*s))/(2*den*s);
   }

   // Leg 2 - Matrix element, summing over two rapidities
   for (unsigned int i=0; i<2; i++) dipg += Born_Parton(M2,-2.*pap1[i],flav1,p)*2./(s-M2-2.*pap3[i]);

   // Output
   return p->alphas()*(1.+s/M2+2.*M2/(s-M2)) * dipg * 96.*pow(M2,2.);
}

// ************************************************************************* //
//  Real quark emission contributions                                        //
// ************************************************************************* //
double Real_Quark_ll(const double s, const double M2, const double pt2, const double cth, const double phi, const unsigned int flav1, const unsigned int flav2, Parameters *p)
{
   return 0.;
}


double Dipole_Quark_ll(const double s, const double M2, const double pt2, const double cth, const double phi, const unsigned int flav1, const unsigned int flav2, Parameters *p)
{
   return 0.;
}


