// ************************************************************************* //
// Dipoles                                                                   //
//                                                                           //
// By Benjamin Fuks - 18.07.2022                                             //
// ************************************************************************* //

// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers ----------------------------------------------- //
#include <cmath>              // Mathematical functions                      //
#include <gsl/gsl_sf_dilog.h> // GSL special functions                       //
// ------------------------------------------------------------------------- //


// ************************************************************************* //
// Regularised splitting functions + residues                                //
// ************************************************************************* //
double Pqq_reg(const double &x) { return -4./3. * (1.+x);                 }
double Pgq_reg(const double &x) { return (pow(x,2.) + pow(1.-x,2.)) / 2.; }

double Phatqq_prime(const double &x) { return 4./3.*(1.-x); }
double Phatgq_prime(const double &x) { return x*(1.-x);     }

// ************************************************************************* //
// Splitting functions                                                       //
// no_delta = withtout the delta(1-x) and plus distribution parts            //
// delta    = delta(1-x) parts + the integration from 0 to xmin of the plus  //
//            distribution component so that the latter is integrated from 0 //
//            to 1)                                                          //
// plus     = plus distribution parts (to be integrated from 0 to 1)         //
// ************************************************************************* //
double Pqq_nodelta(const double &x) { return Pqq_reg(x) + 8./3. / (1.-x); }
double Pqq_delta  (const double &x) { return 2. + 8./3.*log(1.-x);        }
double Pqq_plus   (const double &x) { return 8./3. / (1.-x);              }

double Pgq_nodelta(const double &x) { return Pgq_reg(x); }


// ************************************************************************* //
// Collinear K functions                                                     //
// ************************************************************************* //
double Kbarqq_nodelta(const double &x) { return Pqq_reg(x)*log((1.-x)/x) + Phatqq_prime(x) + 8./3.*log((1.-x)/x)/(1.-x);       }
double Kbarqq_delta  (const double &x) { return -4./3.*(5. - 2./3.*pow(M_PI,2.) - 2.*gsl_sf_dilog(1.-x) - pow(log(1.-x),2.) ); }
double Kbarqq_plus   (const double &x) { return 8./3.*log((1.-x)/x)/(1.-x);                                                    }

double Kbargq_nodelta(const double &x) { return Pgq_reg(x)*log((1.-x)/x) + Phatgq_prime(x); }

double Kqq_nodelta(const double &x) { return Kbarqq_nodelta(x)  + Pqq_reg(x)*log(1.-x) + 8./3.*log(1.-x)/(1.-x); }
double Kqq_delta  (const double &x) { return Kbarqq_delta(x) + 4./3. * (pow(log(1.-x),2.) - pow(M_PI,2.)/3.);    }
double Kqq_plus   (const double &x) { return Kbarqq_plus(x) + 8./3.*log(1.-x)/(1.-x);                            }

double Kgq_nodelta(const double &x) { return Kbargq_nodelta(x)  + Pgq_reg(x)*log(1.-x); }

