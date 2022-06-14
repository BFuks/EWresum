// ************************************************************************* //
// Virtual contributions to the cross section                                //
//                                                                           //
// By Benjamin Fuks - 14.06.2022                                             //
// ************************************************************************* //



// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
#include <iostream>
// -------- Classes -------------------------------------------------------- //
#include "parameters.h"     // All the parameters of the calculation         //
#include "clooptools.h"     // Looptools                                     //
// ------------------------------------------------------------------------- //

static void SetC(double s, int ieps, std::complex<double> *c0)
{
  setlambda(-ieps);
  *c0 = C0i(cc0, 0., 0., s, 0., 0., 0.); 
  exit(1);
}

// ************************************************************************* //
//  Virtual corrections.                                                     //
// ************************************************************************* //
double virt_ll(const double s, const double t, Parameters *p)
{
   // PV integrals
   std::cout << s << " " << t << std::endl;
   int ieps=0;
   std::complex<double> c0;
   SetC(s,ieps,&c0);

   // result (including colour, flux and phase space)
   return 0.;
}

// ************************************************************************* //
//  Integrated dipole contributions.                                         //
// ************************************************************************* //
double dipole_ll(const double s, const double t, Parameters *p)
{
   // result (including colour, flux and phase space)
   return 0.;
}

