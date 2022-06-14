// ************************************************************************* //
// Management of the VEGAS parameters                                        //
//                                                                           //
// By Benjamin Fuks - 12.01.2022                                             //
// ************************************************************************* //


#ifndef VEGAS_CONFIG_H
#define VEGAS_CONFIG_H


// ************************************************************************* //
// Definition of the class Parameters                                        //
// ************************************************************************* //
class VegasConfig
{
  public:
    VegasConfig()  { };
    ~VegasConfig() { };

    // Public accessors
    int VegasVerbose()      { return vegas_verbose_;   }
    double VegasPrecision() { return vegas_precision_; }
    int VegasMaxIter()      { return vegas_max_iter_;  }

    // Mutators
    void SetVegasVerbose(int val)      { vegas_verbose_   = val; }
    void SetVegasMaxIter(int val)      { vegas_max_iter_  = val; }
    void SetVegasPrecision(double val) { vegas_precision_ = val; }


  private:
    // Class members
    int vegas_verbose_, vegas_max_iter_;
    double vegas_precision_;

};
#endif
