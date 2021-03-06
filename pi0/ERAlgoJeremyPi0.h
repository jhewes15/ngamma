/**
 * \file ERAlgoJeremyPi0.h
 *
 * \ingroup Algo
 * 
 * \brief Class def header for a class ERAlgoJeremyPi0
 *
 * @author jhewes15
 */

/** \addtogroup Algo

    @{*/

#ifndef ERTOOL_ERALGOJEREMYPI0_H
#define ERTOOL_ERALGOJEREMYPI0_H

#include "ERTool/Base/AlgoBase.h"
#include "NGammaBase.h"

namespace ertool {

  /**
     \class ERAlgoJeremyPi0
     User custom Algorithm class made by jeremy h
   */
  class ERAlgoJeremyPi0 : public AlgoBase {
  
  public:

    /// Default constructor
    ERAlgoJeremyPi0(const std::string& name="Pi0");

    /// Default destructor
    virtual ~ERAlgoJeremyPi0(){};
    
    /// Set verbosity
    void SetVerbose(bool on){ _verbose = on; _nGamma->SetVerbose(on); };
    
    /// Enable training mode
    void TrainingMode(bool on){ _trainingMode = on; };

    /// Reset function
    void Reset();

    /// Function to accept fclite::PSet
    void AcceptPSet(const ::fcllite::PSet& cfg);

    /// Called @ before processing the first event sample
    void ProcessBegin();

    /// Function to evaluate input showers and determine a score
    bool Reconstruct(const EventData &data, ParticleGraph& graph);

    /// Called after processing the last event sample
    void ProcessEnd(TFile* fout=nullptr);
    
    /// Take identified pi0s and add them to ParticleGraph
    void AddPi0s(ParticleGraph& graph, CombinationScoreSet_t particles);
    
    void SetCut(int n){ _cut = (double)n / 10; };
    
  protected:
    
    std::string _name;

    bool _verbose      = false;
    bool _trainingMode = false;
    
    double _cut = 0.7;
    
    NGammaBase *_nGamma;
  };
}
#endif

/** @} */ // end of doxygen group 
