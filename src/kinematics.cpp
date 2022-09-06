// ************************************************************************* //
// Kinematic computations                                                    //
//                                                                           //
// By Benjamin Fuks - 06.09.2022                                             //
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
// Two-body phase space                                                      //
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


// ************************************************************************* //
// Three-body phase space in the collinear limit (for P+K integration)       //
// ************************************************************************* //
void Kinematics2to2PK(double &xa, double &xb, double &z, double *t, double *jac, double *s, double *x, ProcessConfig *p)
{
   // Computes xa + integration from 0->1 becomes from tau to 1
   xa=pow(p->tauh(),1.-x[0]);

   // Computes xb and z
   if(p->M()==-1)
   {
     xb=pow(p->tauh()/xa,1.-x[3]);
     z = pow(p->tauh()/(xa*xb),1.-x[1]);
   }
   else
   {
     z = pow(p->tauh()/xa,1.-x[1]);
     xb=p->tauh()/(z*xa);
   }

   // Computes { shat, shatilde }
   s[0] = xa*xb*p->sh();
   s[1] = xa*xb*z*p->sh();

   // Computes t/ttilde + integration from 0->1 becomes from tmin to tmax 
   double sqkal[2] = { sqrt(Kallen(s[0],p->m1sq(),p->m2sq())), sqrt(Kallen(s[1],p->m1sq(),p->m2sq())) };
   for (unsigned int i=0; i<2;i++)
   {
      double tmin = (-s[i] + p->m2sq() + p->m1sq() - sqkal[i])/2.;
      double tmax = tmin + sqkal[i];
      t[i]=(tmax-tmin)*x[2]+tmin;
   }

   // Jacobian divided by (xa xb) cf. PDF (includes xb*z or xb factor for dsig/dM)
   // 0 = colilinear phase space
   // 1 = with the + distribution
   // 2 =  for the delta(1-z)
   if(p->M()==-1)
   {
     jac[0] = -sqkal[1]*z*log(p->tauh())*log(p->tauh()/xa)*log(p->tauh()/(xa*xb));
     jac[1] = -sqkal[0]*z*log(p->tauh())*log(p->tauh()/xa)*log(p->tauh()/(xa*xb));
     jac[2] =  sqkal[0]*log(p->tauh())*log(p->tauh()/xa);
   }
   else
   {
     jac[0] =  sqkal[1]*z*log(p->tauh())*log(p->tauh()/(xa*xb));
     jac[1] =  sqkal[1]*z*log(p->tauh())*log(p->tauh()/(xa*xb));
     jac[2] = -sqkal[1]*z*log(p->tauh());
   }
}



// ************************************************************************* //
// Three-body phase space                                                    //
// ************************************************************************* //
void Kinematics2to3(double &xa, double &xb, double &M2, double &pt2, double &cthp, double &php, double &jac, double &s, double *x, ProcessConfig* p)
{
   // Computes xa + integration from 0->1 becomes from tau to 1
   xa=pow(p->tauh(),1.-x[0]);

   // Computes xb
   xb=pow(p->tauh()/xa,1.-x[1]);

   // Computes M2 and shat
   s = xa*xb*p->sh();
   if(p->M()==-1) M2 = p->M2min() * pow(s/p->M2min(),x[5]);
   else           M2=p->M2();

   // Computes pt2
   double pt2max = pow(s-M2, 2.) / (4.*s);
   pt2 = pt2max*x[2];

   // Computes cos(theta') and phi'
   cthp = 2.*x[3]-1.;
//   double thp = M_PI * x[3];
//   cthp = cos(thp);
   php = M_PI * x[4];

   // Jacobian divided by (xa xb) cf. PDF (includes xb factor for dsig/dM)
   //  + factor of 2 for the phi' integral that is only done half
   if(p->M()==-1) jac = 4.*M_PI * M2 * pt2max * log(p->tauh()) * log(p->tauh()/xa) * log(s/p->M2min());
//   if(p->M()==-1) jac = 2.*M_PI * M_PI * M2 * pt2max * log(p->tauh()) * log(p->tauh()/xa) * log(s/p->M2min()) * sin(thp);
   else           jac = 4.*M_PI * pt2max * log(p->tauh()) * log(p->tauh()/xa);

   return;

}


