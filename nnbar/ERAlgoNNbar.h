/**
 * \file ERAlgoNNbar.h
 *
 * \ingroup nnbar
 * 
 * \brief Class def header for a class ERAlgoNNbar
 *
 * @author jeremy h
 */

/** \addtogroup nnbar

    @{*/

#ifndef ERTOOL_ERALGONNBAR_H
#define ERTOOL_ERALGONNBAR_H

#include "ERTool/Base/AlgoBase.h"
#include "NGammaBase.h"

namespace ertool {

  /**
     \class ERAlgoNNbar
     User custom Algorithm class made by jeremy h
   */
  class ERAlgoNNbar : public AlgoBase {
  
  public:

    /// Default constructor
    ERAlgoNNbar(const std::string& name="NNbar");

    /// Default destructor
    virtual ~ERAlgoNNbar(){};
    
    /// Set verbosity
    void SetVerbose(bool on){ _verbose = on; _nGamma->SetVerbose(on); };
    
    /// Set training mode (true = enabled)
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
