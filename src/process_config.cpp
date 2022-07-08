// ************************************************************************* //
// Management of the process configuration                                   //
//                                                                           //
// By Benjamin Fuks - 14.06.2022                                             //
// ************************************************************************* //

// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers ----------------------------------------------- //
#include <cmath>            // Mathematical functions                        //
#include <complex>          // Complex numbers library                       //
using namespace std::complex_literals;
// -------- Classes -------------------------------------------------------- //
#include "process_config.h"     // Parameters                                //
// -------- Functions ------------------------------------------------------ //
//                                                                           //
// ------------------------------------------------------------------------- //


// ************************************************************************* //
//  Checking the parameter consistency                                       //
// ************************************************************************* //
void C_init(double C[6][6])
{
  for(unsigned int i=0; i<6; i++)
    for(unsigned int j=0; j<6; j++)
      C[i][j]=4.;
}

bool ProcessConfig::Check()
{
  // Printing information and initialisastion
  info("Checking the consistency of the process parameters");

  // Masses
  switch(std::abs(PDG1()))
  {
    case 11: flavour1_ = 0; m1_=.0; m1sq_=0.; break;
    case 12: flavour1_ = 1; m1_=.0; m1sq_=0.; break;
    case 13: flavour1_ = 2; m1_=.0; m1sq_=0.; break;
    case 14: flavour1_ = 3; m1_=.0; m1sq_=0.; break;
    case 15: flavour1_ = 4; m1_=.0; m1sq_=0.; break;
    case 16: flavour1_ = 5; m1_=.0; m1sq_=0.; break;
    default: m1_=.0; m1sq_=0.;
  }

  switch(std::abs(PDG2()))
  {
    case 11: flavour2_ = 0; m2_=.0; m2sq_=0.; break;
    case 12: flavour2_ = 1; m2_=.0; m2sq_=0.; break;
    case 13: flavour2_ = 2; m2_=.0; m2sq_=0.; break;
    case 14: flavour2_ = 3; m2_=.0; m2sq_=0.; break;
    case 15: flavour2_ = 4; m2_=.0; m2sq_=0.; break;
    case 16: flavour2_ = 5; m2_=.0; m2sq_=0.; break;
    default: m2_=.0; m2sq_=0.;
  }

  // Hadronic center of mass energy
  sh_   = pow(EBeam1()+EBeam2(),2.);
  if(M_==-1) tauh_ = pow(std::max(m1_+m2_, Mmin_),2.)/sh_;
  else       tauh_ = pow(M_,2.)/sh_;

  // Couplings
  // EW inputs
  ee2_   = 4.*M_PI/aewm1_;  ee_ = sqrt(ee2_);

  mzsq_  = pow(mz_,2);
  cmzsq_ = mzsq_ + 1i*mz_*gamz_;
  mwsq_  = mzsq_/2. + sqrt(pow(mz_,4.)/4. - M_PI*mzsq_/(sqrt(2.)*gf_*aewm1_));
  mw_   = sqrt(mwsq_);
  cmwsq_ = mwsq_ + 1i*mw_*gamw_;

  sw2_  = 1.-mwsq_/mzsq_;
  cw_   = sqrt(1.-sw2_);
  sw_   = sqrt(sw2_);

  // d u s c b t (LHAPDF ordering)
  C_init(Aqq_); C_init(All_);
  C_init(ZLqq_); C_init(ZRqq_); C_init(ZLll_); C_init(ZRll_);
  Aqq_[0][0] = ee_/3.;  Aqq_[1][1] = -2.*ee_/3.;
  Aqq_[2][2] = ee_/3.;  Aqq_[3][3] = -2.*ee_/3.;
  Aqq_[4][4] = ee_/3.;  Aqq_[5][5] = -2.*ee_/3.;

  All_[0][0] = ee_;     All_[1][1] = 0.;
  All_[2][2] = ee_;     All_[3][3] = 0.;
  All_[4][4] = ee_;     All_[5][5] = 0.;

  ZLqq_[0][0] = ee_/(sw_*cw_)*(sw2_/3.-.5);  ZLqq_[1][1] = ee_/(sw_*cw_)*(-2.*sw2_/3.+.5);
  ZLqq_[2][2] = ee_/(sw_*cw_)*(sw2_/3.-.5);  ZLqq_[3][3] = ee_/(sw_*cw_)*(-2.*sw2_/3.+.5);
  ZLqq_[4][4] = ee_/(sw_*cw_)*(sw2_/3.-.5);  ZLqq_[5][5] = ee_/(sw_*cw_)*(-2.*sw2_/3.+.5);
  ZRqq_[0][0] = ee_/(sw_*cw_)*(sw2_/3.   );  ZRqq_[1][1] = ee_/(sw_*cw_)*(-2.*sw2_/3.   );
  ZRqq_[2][2] = ee_/(sw_*cw_)*(sw2_/3.   );  ZRqq_[3][3] = ee_/(sw_*cw_)*(-2.*sw2_/3.   );
  ZRqq_[4][4] = ee_/(sw_*cw_)*(sw2_/3.   );  ZRqq_[5][5] = ee_/(sw_*cw_)*(-2.*sw2_/3.   );

  ZLll_[0][0] = ee_/(sw_*cw_)*(sw2_   -.5);  ZLll_[1][1] = ee_/(sw_*cw_)*(           +.5);
  ZLll_[2][2] = ee_/(sw_*cw_)*(sw2_   -.5);  ZLll_[3][3] = ee_/(sw_*cw_)*(           +.5);
  ZLll_[4][4] = ee_/(sw_*cw_)*(sw2_   -.5);  ZLll_[5][5] = ee_/(sw_*cw_)*(           +.5);
  ZRll_[0][0] = ee_/(sw_*cw_)*(sw2_      );  ZRll_[1][1] = 0.;
  ZRll_[2][2] = ee_/(sw_*cw_)*(sw2_      );  ZRll_[3][3] = 0.;
  ZRll_[4][4] = ee_/(sw_*cw_)*(sw2_      );  ZRll_[5][5] = 0.;

  // Exit
  return true;
}
