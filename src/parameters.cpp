// ************************************************************************* //
// Management of the calculation parameters                                  //
//                                                                           //
// By Benjamin Fuks - 08.07.2022                                             //
// ************************************************************************* //

// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers ----------------------------------------------- //
#include <fstream>        // Read/write files                                //
#include <sstream>        // String streams                                  //
#include <string>         // Strings                                         //
#include <vector>         // Vectors                                         //
// -------- Classes -------------------------------------------------------- //
#include "messages.h"       // Message services                              //
#include "parameters.h"     // All the parameters of the calculation         //
// -------- Functions ------------------------------------------------------ //
std::vector<std::string> tokenise(const std::string&, const std::string&);   //
// ------------------------------------------------------------------------- //

// ************************************************************************* //
//  Constructor                                                              //
// ************************************************************************* //
Parameters::Parameters(bool &dbg,  const std::string &s1)
{
  // Printing information and initialisastion
  info("Initialising the parameters of the calculation");
  vegas_   = new VegasConfig();
  process_ = new ProcessConfig();

  // Debug mode
  SetDebug(dbg);

  // Checking whether file exists
  std::ifstream myfile(s1.c_str());
  if(!myfile.is_open()) error("Invalid configuration file");

  // Decoding information
  std::string line;
  while(myfile.good())
  {
    std::getline(myfile,line);

    // Safety
    if(line.length()==0) continue;

    // Comments
    const std::string whitespace = " \n\r\t\f\v";
    size_t start = line.find_first_not_of(whitespace);
    if(line.substr(start).rfind("#", 0)==0) continue;

    // We have a non-empty line  -> getting key and values
    std::vector<std::string> decoded_line = tokenise(line.substr(start),"#");
    decoded_line = tokenise(decoded_line[0]," ");
    if(decoded_line.size()!=2) error("Error in the config file format (parameter value # comment): \n              " + line);
    std::istringstream valuestream(decoded_line[1]);

    // Calculation configuration
    if     (decoded_line[0]=="vegas_verbosity") { int    ival; valuestream >> ival; SetVegasVerbose(ival);   }
    else if(decoded_line[0]=="vegas_max_iter")  { int    ival; valuestream >> ival; SetVegasMaxIter(ival);   }
    else if(decoded_line[0]=="vegas_precision") { double dval; valuestream >> dval; SetVegasPrecision(dval); }
    else if(decoded_line[0]=="lhapdf_set_born") { unsigned int ival; valuestream >> ival; SetLHAPDFSetBorn(ival); }
    else if(decoded_line[0]=="lhapdf_set")      { unsigned int ival; valuestream >> ival; SetLHAPDFSet(ival); }
    else if(decoded_line[0]=="lhapdf_verbosity"){ int    ival; valuestream >> ival; SetLHAPDFVerbose(ival); }

    // Process
    else if(decoded_line[0]=="beam1_energy"){ double dval; valuestream >> dval; SetEBeam1(dval); }
    else if(decoded_line[0]=="beam2_energy"){ double dval; valuestream >> dval; SetEBeam2(dval); }
    else if(decoded_line[0]=="final_pdgid_1"){ int   ival; valuestream >> ival; SetPDG1(ival); }
    else if(decoded_line[0]=="final_pdgid_2"){ int   ival; valuestream >> ival; SetPDG2(ival); }
    else if(decoded_line[0]=="M")      { double dval; valuestream >> dval; SetM(dval);    }
    else if(decoded_line[0]=="Mmin")   { double dval; valuestream >> dval; SetMmin(dval); }
    else if(decoded_line[0]=="muR")    { double dval; valuestream >> dval; SetmuR(dval);  }
    else if(decoded_line[0]=="muF")    { double dval; valuestream >> dval; SetmuF(dval);  }

    // parameters
    else if(decoded_line[0]=="aewm1"){ double dval; valuestream >> dval; SetaEWM1(dval); }
    else if(decoded_line[0]=="gf")   { double dval; valuestream >> dval; SetGF(dval); }
    else if(decoded_line[0]=="mz")   { double dval; valuestream >> dval; SetMZ(dval); }
    else if(decoded_line[0]=="gz")   { double dval; valuestream >> dval; SetGZ(dval); }
    else if(decoded_line[0]=="gw")   { double dval; valuestream >> dval; SetGW(dval); }

    else error("Unknown keyword in the config file: " + decoded_line[0]);

    if(Debug()) debug("Parameter " + decoded_line[0] + " = " + decoded_line[1]);
  }

  // Done
  myfile.close();
  info("  --> done!");
}



// ************************************************************************* //
//  Checking the parameter consistency                                       //
// ************************************************************************* //
bool Parameters::Check()
{
  // Printing information and initialisastion
  info("Checking the consistency of the input file");

  // Checking whether the selected parton density set is avaialble
  LHAPDF::pair<std::string, int> set_id = LHAPDF::lookupPDF(PDFId());
  LHAPDF::pair<std::string, int> set_idb= LHAPDF::lookupPDF(PDFIdBorn());
  bool available_pdf = false, available_pdfb = false;
  for (const std::string& s : LHAPDF::availablePDFSets())
  {
    if(s.compare(set_id.first)==0)
    {
      if(Debug()) debug("The PDF set " + set_id.first + " has been found");
      available_pdf = true;
    }
    if(s.compare(set_idb.first)==0)
    {
      if(Debug()) debug("The PDF set " + set_idb.first + " has been found");
      available_pdfb = true;
    }
  }
  if(!available_pdf)  error("The PDF set " + set_id.first + " is not available");
  if(!available_pdfb) error("The PDF set " + set_idb.first + " is not available");

  // Collision parameters
  Process()->Check();


  // Everything is normal
  info("  --> done!");
  return true;
}
