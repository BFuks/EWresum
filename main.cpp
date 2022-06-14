// ************************************************************************* //
// Main code for resummed cross section calculations                         //
//                                                                           //
// By Benjamin Fuks - 12.01.2022                                             //
// ************************************************************************* //



// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Classes ---------------------------------------------- //
#include "messages.h"     // Message services                      //
#include "parameters.h"   // All the parameters of the calculation //
// -------- Functions -------------------------------------------- //
void DecodeArgs(int, char**, bool&); // Arguments from the shell   //
void DisplayXsec(const double&, const double&, const std::string&);//
void Integrate(double (*)(double *,size_t,void *),double&,double&,size_t,Parameters*);//
double Born(double*,size_t,void*);                                 //
double NLO(double*,size_t,void*);                                  //
// --------------------------------------------------------------- //



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

  // Main integrator
  info("Born cross section calculation");
  Params->LHAPDF_init(Params->PDFIdBorn());
  double res0=0., err0=0.;
  size_t ndims = 2;
  if(Params->M()==-1) { ndims+=1; }
  Integrate(&Born,res0,err0,ndims,Params);

  // Display of the results
  info("Results of the calculation");
  DisplayXsec(res0, err0, "  --> Final");

  // End of the program
  return 0;
}
