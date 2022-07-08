// ************************************************************************* //
// Management of the process configuration                                   //
//                                                                           //
// By Benjamin Fuks - 08.07.2022                                             //
// ************************************************************************* //


#ifndef PROCESS_CONFIG_H
#define PROCESS_CONFIG_H
// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers ----------------------------------------------- //
#include <complex>          // Complex numbers library                       //
// -------- Classes -------------------------------------------------------- //
#include "messages.h"       // Message services                              //
// ------------------------------------------------------------------------- //



class ProcessConfig
{
  public:
    // Constructors
    ProcessConfig() { };
    ~ProcessConfig() { };

    // Checking the parameter consistency
    bool Check();

    // Public accessors
    double EBeam1()      { return ebeam1_;   }
    double EBeam2()      { return ebeam2_;   }
    int PDG1()           { return pdgid1_;   }
    int PDG2()           { return pdgid2_;   }
    unsigned int flav1() { return flavour1_; }
    unsigned int flav2() { return flavour2_; }

    double tauh()   { return tauh_;   }
    double sh()     { return sh_;     }
    double M()      { return M_;      }
    double m1sq()   { return m1sq_;   }
    double m2sq()   { return m2sq_;   }

    double muR()    { return muR_;    }
    double muR2()   { return muR2_;   }
    double muF()    { return muF_;    }
    double muF2()   { return muF2_;   }

    double mz()     { return mz_;     }
    double mz2()    { return mzsq_;   }
    std::complex<double> cmz2() { return cmzsq_;  }

    double ee()     { return ee_;     }
    double ee2()    { return ee2_;    }
    double  Aqq(unsigned int i, unsigned int j) {return Aqq_[i][j];}
    double  All(unsigned int i, unsigned int j) {return All_[i][j];}
    double ZLqq(unsigned int i, unsigned int j) {return ZLqq_[i][j];}
    double ZRqq(unsigned int i, unsigned int j) {return ZRqq_[i][j];}
    double ZLll(unsigned int i, unsigned int j) {return ZLll_[i][j];}
    double ZRll(unsigned int i, unsigned int j) {return ZRll_[i][j];}

    // Public Mutators
    void SetEBeam1(double val)   { ebeam1_ = val; }
    void SetEBeam2(double val)   { ebeam2_ = val; }
    void SetPDG1(int val)        { pdgid1_ = val; }
    void SetPDG2(int val)        { pdgid2_ = val; }
    void SetaEWM1(double val)    { aewm1_  = val; }
    void SetGF(double val)       { gf_     = val; }
    void SetMZ(double val)       { mz_     = val; }
    void SetGZ(double val)       { gamz_  = val; }
    void SetGW(double val)       { gamw_   = val; }
    void SetM(double val)        { M_      = val; }
    void SetMmin(double val)     { Mmin_   = val; }
    void SetmuR(double val)      { muR_    = val; muR2_ = pow(val, 2.); }
    void SetmuF(double val)      { muF_    = val; muF2_ = pow(val, 2.); }

  private:
    // class members
    double ebeam1_, ebeam2_;  // Beam energies
    int pdgid1_, pdgid2_;     // Final-state particle PDG codes
    unsigned int flavour1_;   // Final-state flavour flag 1
    unsigned int flavour2_;   // Final-state flavour flag 2

    // Internal parameters
    double sh_;               // Hadronic centre of mass energy
    double M_, Mmin_;         // Invariant mass
    double tauh_;             // Reduced invariant mass
    double muR_, muR2_;       // Renormalisation scale
    double muF_, muF2_;       // Factorisation scale

    // masses
    double mz_, mzsq_, mw_, mwsq_;       // Weak boson masses
    double gamz_, gamw_;                 // Weak boson widths
    std::complex<double> cmzsq_, cmwsq_; // Weak boson complex masses
    double m1_, m2_, m1sq_, m2sq_;       // Final state masses

    // Couplings
    double aewm1_, gf_;                         // EW inputs
    double ee_, ee2_, cw_, sw2_, sw_;           // EW parameters
    double Aqq_[6][6], ZLqq_[6][6], ZRqq_[6][6];  // V-qq couplings
    double All_[6][6], ZLll_[6][6], ZRll_[6][6];  // V-ll couplings

};
#endif
