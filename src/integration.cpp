// ************************************************************************* //
// Integration main routines                                                  //
//                                                                           //
// By Benjamin Fuks - 14.06.2022                                             //
// ************************************************************************* //



// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers ------------------------------------- //
#include <sstream>        // String streams                        //
// -------- Classes ---------------------------------------------- //
#include "messages.h"     // Message services                      //
#include "parameters.h"   // All the parameters of the calculation //
#include <gsl/gsl_monte_vegas.h>  // GSL                           //
// -------- Functions -------------------------------------------- //
void DisplayXsec(const double&, const double&, const std::string&);//
// --------------------------------------------------------------- //


// ************************************************************************* //
//  Main integration routine, steering vegas.                                //
// ************************************************************************* //
void Integrate(double (*fctn)(double *x, size_t dim, void *jj), double& res, double& err, size_t ndims, Parameters* Params)
{
  // Initializing the random number generator 
  info("  --> Initialisation of the integrator");
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  time_t seed = time(NULL);
  gsl_rng_set(r,seed);

  // Integrand
  gsl_monte_function I; I.f = fctn; I.dim=ndims; I.params=Params;

  // Number of integrations and calls
  size_t calls=10000;

  // Integration bounds
  double xmin[I.dim], xmax[I.dim];
  for(size_t i=0; i<I.dim; i++) { xmin[i]=0.; xmax[i]=1.; }

  // Integration
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(I.dim);

  // Warm-up
  s->stage=0;  s->iterations=5; s->verbose=Params->Vegas()->VegasVerbose();
  gsl_monte_vegas_integrate(&I, xmin, xmax, I.dim, calls/50, r, s, &res, &err);
  DisplayXsec(res,err,"  --> Warm-up"); 

  // Real things: stops if the precision reaches 0.1% or if one oscillates
  double prec=1e9; int counter=1; s->stage=1;
  while(prec>Params->Vegas()->VegasPrecision() && counter<=Params->Vegas()->VegasMaxIter())
  {
    // Integral
    gsl_monte_vegas_integrate(&I, xmin, xmax, I.dim, calls, r, s, &res, &err);

    // Display result
    std::ostringstream ocnt; ocnt << counter; std::string cntstr= ocnt.str();
    DisplayXsec(res,err,"  --> Refine-"+cntstr);

    // Preparing next iteration
    prec=std::abs(err/res); counter++; calls*=5; 
  }

  // Cleaning the memory and closing the file
  gsl_monte_vegas_free(s);
  gsl_rng_free(r);
}

