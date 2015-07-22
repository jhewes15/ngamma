/**
 * \file NGammaBase.h
 *
 * \ingroup Algo
 * 
 * \brief Class def header for a class NGammaBase
 *
 * @author jeremy h
 */

/** \addtogroup Algo

    @{*/
#ifndef ERNGAMMABASE_H
#define ERNGAMMABASE_H

#include "ERTool/Algo/AlgoEMPart.h"

/**
   \class NGammaBase
   User defined class NGammaBase ... these comments are used to generate
   doxygen documentation!
 */

namespace ertool {

  typedef std::pair<Combination_t,double> CombinationScore_t;
  typedef std::vector<CombinationScore_t> CombinationScoreSet_t;

  class NGammaBase{

  public:

    /// Default constructor
    NGammaBase();

    /// Default destructor
    virtual ~NGammaBase(){};
    
    /// function to control verbosity
    void SetVerbose(bool on) { _verbose = on; }
    
    ::fcllite::PSet& OutputPSet() { return _pset; }
    
    /// Function to accept fclite::PSet
    void AcceptPSet(const ::fcllite::PSet& cfg);
    
    /// Loop over shower combinations and set score
    CombinationScoreSet_t GammaComparison(const ParticleGraph& graph, int n = 0);
    
    /// Find best combinations in event
    CombinationScoreSet_t EventSelection(CombinationScoreSet_t candidates, double threshold);
    
    /// Likelihood function to score a combination of gammas
    double Likelihood(ParticleGraph graph, Combination_t combination);
    
    /// Functions to train PDFs
    void FillEnergyPdf(double energy);
    void FillMomentumPdf(geoalgo::Vector momentum);
    void FillMassPdf(double mass);
    fcllite::PSet GetParams();
    
    /// Function to set a RooFit parameter
    void SetParam(RooRealVar* var, double val);
    
    /// Functions to set external variables
    void SetExternalEnergy(double energy){ _extEnergy = energy; };
    void SetExternalMomentum(geoalgo::Vector momentum){ delete _extMomentum; _extMomentum = new geoalgo::Vector(momentum); };
    
  protected:
    
    bool _verbose;
    
    // parameter set
    ::fcllite::PSet _pset;
    
    // pdf things
    PdfFactory  _factory;
    
    RooRealVar* _energyVar;
    RooAbsPdf*  _energyPdf;
    RooDataSet* _energyData;
    
    RooRealVar* _momentumVar;
    RooAbsPdf*  _momentumPdf;
    RooDataSet* _momentumData;
    
    RooRealVar* _massVar;
    RooAbsPdf*  _massPdf;
    RooDataSet* _massData;
    
    // variables for pdfs
    bool _useEnergyPdf   = false;
    bool _useMomentumPdf = false;
    bool _useMassPdf     = false;
    
    // external variables for likelihood
    double _extEnergy;
    geoalgo::Vector *_extMomentum;
  };
}
#endif
/** @} */ // end of doxygen group 

