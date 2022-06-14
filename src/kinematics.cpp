// ************************************************************************* //
// Kinematic computations                                                    //
//                                                                           //
// By Benjamin Fuks - 14.06.2022                                             //
// ************************************************************************* //

// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers ----------------------------------------------- //
#include <cmath>            // Mathematical functions                        //
// -------- Classes -------------------------------------------------------- //
#include "process_config.h" // Class containing the process details          //
// ------------------------------------------------------------------------- //



// ************************************************************************* //
// Kallen function                                                           //
// ************************************************************************* //
double Kallen(double x, double y, double z)
  { return ( pow(x,2.)+pow(y,2.)+pow(z,2.)-2.*(x*y+y*z+x*z) ); }


// ************************************************************************* //
// Kinematics                                                                //
// ************************************************************************* //
void Kinematics2to2(double &xa, double &xb, double &t, double &jac, double &s, double *x, ProcessConfig *p)
{
   // Computes xa + integration from 0->1 becomes from tau to 1
   xa=pow(p->tauh(),1.-x[0]);

   // Computes xb
   if(p->M()==-1) xb=pow(p->tauh()/xa,1.-x[2]);
   else           xb=p->tauh()/xa;

   // Computes shat
   s = xa*xb*p->sh();

   // Computes t + integration from 0->1 becomes from tmin to tmax 
   double sqkal = sqrt(Kallen(s,p->m1sq(),p->m2sq()));
   double tmin = (-s + p->m2sq() + p->m1sq() - sqkal)/2.;
   double tmax = tmin + sqkal;
   t=(tmax-tmin)*x[1]+tmin;

   // Jacobian divided by (xa xb) cf. PDF (includes xb factor for dsig/dM)
   if(p->M()==-1) jac =  sqkal*log(p->tauh())*log(p->tauh()/xa);
   else           jac = -sqkal*log(p->tauh());
}

