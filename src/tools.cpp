// ************************************************************************* //
// Various tools                                                             //
//                                                                           //
// By Benjamin Fuks - 14.06.2022                                             //
// ************************************************************************* //


// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers ------------------------------------- //
#include <cmath>            // Mathematical functions           //
#include <sstream>    // String streams                            //
#include <string>     // Strings                                   //
#include <vector>     // Vectors                                   //
// -------- Classes ---------------------------------------------- //

// -------- Classes ---------------------------------------------- //
#include "messages.h"     // Message services                      //
// --------------------------------------------------------------- //



// ************************************************************************* //
//  Decoing the arguments passed from the shell                              //
// ************************************************************************* //
void DecodeArgs(int argc, char* argv[], bool &debug)
{
  // Number of mandatory arguments
  unsigned int nargs = 0;

  // Loop over all passed arguments
  for (unsigned int i=1; i<argc; i++)
  {
    std::string str = argv[i];
    // mandatory argument
    if(str.rfind("--", 0)!=0 && str.rfind("-", 0)!=0) nargs++;
    // options
    else if (str.compare("--debug")==0 || str.compare("-d")==0 || str.compare("-D")==0) debug = true;
    else
    {
      std::string msg("Unknown option: "), msg2(argv[i]);
      error(msg.append(msg2).c_str());
    }
  }

  // Safety check
  if(nargs!=1) error("This function requires 1 argument: the path to a configuration file");
}



// ************************************************************************* //
// Printing a xsection                                                       //
// ************************************************************************* //
void DisplayXsec(const double &r, const double &e, const std::string &tag)
{
  // From double to strings
  std::ostringstream osr; osr << r*0.38937966e9; std::string rstr= osr.str();
  std::ostringstream ose; ose << e*0.38937966e9; std::string estr= ose.str();

  // Precision of the result
  double prec=0.;
  if(r!=0.) prec=std::abs(100.*e/r);
  std::ostringstream opr; opr.precision(2); 
  opr << std::fixed << prec; std::string prstr= opr.str();

  // Printing
  info(tag + " results: " + rstr + " +/- " + estr + " pb  (@ " + prstr + "%)");
}


// ************************************************************************* //
// Printing a xsection                                                       //
// ************************************************************************* //
std::vector<std::string> tokenise(const std::string &str, const std::string &token=" ")
{
    // Intialisation
    int start = 0, end = str.find(token);
    std::vector<std::string>result;

    // Tokenisation
    while (end != -1)
    {
      if(str.substr(start, end - start).size()>0) result.push_back(str.substr(start, end - start));
      start = end + token.size();
      end   = str.find(token, start);
    }
    if(str.substr(start, end - start).size()>0) result.push_back(str.substr(start, end - start));

    // Exit
    return result;
}
