// ************************************************************************* //
// Main code for resummed cross section calculations                         //
//                                                                           //
// By Benjamin Fuks - 12.01.2022                                             //
// ************************************************************************* //



// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Classes -------------------------------------------------------- //
#include "messages.h"     // Message services                                //
#include "parameters.h"   // All the parameters of the calculation           //
#include "clooptools.h"   // Looptools                                       //
// -------- Functions --------------------------------------------           //
void DecodeArgs(int, char**, bool&); // Arguments from the shell             //
void DisplayXsec(const double&, const double&, const std::string&);          //
void Integrate(double (*)(double *,size_t,void *),double&,double&,size_t,Parameters*);//
double Born(double*,size_t,void*);                                           //
double Virt(double*,size_t,void*);                                           //
// ------------------------------------------------------------------------- //



// ************************************************************************* //
//  Main code                                                                //
// ************************************************************************* //
int main(int argc, char* argv[])
{
  // Initializing message services
  InitMessages();

  // Decoding the provided arguments
  bool debug = false;
  DecodeArgs(argc,argv,debug);

  // Loading the model parameters, the collision setup, etc.
  Parameters *Params = new Parameters(debug, argv[1]);
  if(!Params->Check()) error("Parameters inconsistent; fix the input file");

  // LO integrator
  info("Born cross section calculation");
  Params->LHAPDF_init(Params->PDFIdBorn());
  double res=0., err=0.;
  size_t ndims = 2;
  if(Params->M()==-1) { ndims+=1; }
  Integrate(&Born,res,err,ndims,Params);
  DisplayXsec(res, err, "  --> Final");

  // inits looptools
  ltini();

  // NLO integrator
  info("NLO cross section calculation");
  Params->LHAPDF_init(Params->PDFId());
  double r_nlo=0., e_nlo=0.;

  info("  --> LO component");
  Integrate(&Born,res,err,ndims,Params);
  DisplayXsec(res, err, "  --> LO component");
  r_nlo=res; e_nlo=err;

  info("  --> \"Virtual + dipole\" component");
  Integrate(&Virt,res,err,ndims,Params);
  DisplayXsec(res, err, "  --> \"Virtual + dipole\" component");
  r_nlo+=res; e_nlo+=err;

  // End of the program (and clearing looptools)
  ltexi();
  return 0;
}
