// ************************************************************************* //
// Management of the calculation parameters                                  //
//                                                                           //
// By Benjamin Fuks - 14.06.2022                                             //
// ************************************************************************* //


#ifndef PARAMETERS_H
#define PARAMETERS_H
// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers ----------------------------------------------- //
#include <complex>          // Complex numbers library                       //
#include <string>           // Strings                                       //
// -------- Classes -------------------------------------------------------- //
#include "vegas_config.h"   // VEGAS configuration                           //
#include "LHAPDF/LHAPDF.h"  // LHAPDF libraries                              //
#include "process_config.h" // Class containing the process details          //
// ------------------------------------------------------------------------- //

// ************************************************************************* //
// Definition of the class Parameters                                        //
// ************************************************************************* //
class Parameters
{
  public:
    // Constructors
    Parameters()  { vegas_ = new VegasConfig(); process_ = new ProcessConfig(); };
    Parameters(bool&, const std::string&);
    ~Parameters() { };

    // Checking the parameter consistency
    bool Check();

    // Public accessors
    bool Debug()               { return debug_;    }
    unsigned int PDFId()       { return pdfID_;    }
    unsigned int PDFIdBorn()   { return pdfIDBorn_;}
    VegasConfig* Vegas()       { return vegas_;    }
    ProcessConfig* Process()   { return process_;  }
    LHAPDF::PDF* pdf()         { return pdf_;      }

    unsigned int flav1()  { return Process()->flav1(); }
    unsigned int flav2()  { return Process()->flav2(); }

    double M()    { return Process()->M();   }
    double mz()   { return Process()->mz();  }
    double mz2()  { return Process()->mz2(); }
    std::complex<double> cmz2() { return Process()->cmz2(); }

    double ee()   { return Process()->ee();  }
    double ee2()  { return Process()->ee2(); }
    double  Aqq(unsigned int i, unsigned int j) { return Process()->Aqq(i,j);  }
    double  All(unsigned int i, unsigned int j) { return Process()->All(i,j);  }
    double ZLqq(unsigned int i, unsigned int j) { return Process()->ZLqq(i,j); }
    double ZRqq(unsigned int i, unsigned int j) { return Process()->ZRqq(i,j); }
    double ZLll(unsigned int i, unsigned int j) { return Process()->ZLll(i,j); }
    double ZRll(unsigned int i, unsigned int j) { return Process()->ZRll(i,j); }

    // LHAPDF initialisation
    void LHAPDF_init(int label) { pdf_ = LHAPDF::mkPDF(label); }

  private:
    // class members
    bool debug_ ;          // debug flag

    VegasConfig* vegas_;     // VEGAS configuration
    ProcessConfig* process_;   // VEGAS configuration

    unsigned int pdfID_, pdfIDBorn_; // LHAPDF label
    LHAPDF::PDF* pdf_;               // LHAPDF handler


    // Mutators
    void SetDebug(bool dbg)             { debug_=dbg;                      }

    void SetVegasVerbose(int val)       { Vegas()->SetVegasVerbose(val);   }
    void SetVegasMaxIter(int val)       { Vegas()->SetVegasMaxIter(val);   }
    void SetVegasPrecision(double val)  { Vegas()->SetVegasPrecision(val); }

    void SetLHAPDFSet(unsigned int val)     { pdfID_ = val;                }
    void SetLHAPDFSetBorn(unsigned int val) { pdfIDBorn_ = val;            }
    void SetLHAPDFVerbose(int val)      { LHAPDF::getConfig().set_entry("Verbosity", val); }

    void SetEBeam1(double val)          { Process()->SetEBeam1(val); }
    void SetEBeam2(double val)          { Process()->SetEBeam2(val); }
    void SetPDG1(double val)            { Process()->SetPDG1(val);   }
    void SetPDG2(double val)            { Process()->SetPDG2(val);   }

    void SetM(double val)               { Process()->SetM(val);      }
    void SetMmin(double val)            { Process()->SetMmin(val);   }
    void SetaEWM1(double val)           { Process()->SetaEWM1(val);  }
    void SetGF(double val)              { Process()->SetGF(val);     }
    void SetMZ(double val)              { Process()->SetMZ(val);     }
    void SetGZ(double val)              { Process()->SetGZ(val);     }
    void SetGW(double val)              { Process()->SetGW(val);     }

};
#endif
