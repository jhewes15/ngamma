/**
 * \file ERAnaNNbar.h
 *
 * \ingroup nnbar
 * 
 * \brief Class def header for a class ERAnaNNbar
 *
 * @author jeremy h
 */

/** \addtogroup nnbar

    @{*/

#ifndef ERTOOL_ERANANNBAR_H
#define ERTOOL_ERANANNBAR_H

#include "ERTool/Base/AnaBase.h"

namespace ertool {

  /**
     \class ERAnaNNbar
     User custom Analysis class made by kazuhiro
   */
  class ERAnaNNbar : public AnaBase {
  
  public:

    /// Default constructor
    ERAnaNNbar(const std::string& name="NNbar");

    /// Default destructor
    virtual ~ERAnaNNbar(){}

    /// Reset function
    virtual void Reset();

    /// Function to accept fclite::PSet
    void AcceptPSet(const ::fcllite::PSet& cfg);

    /// Called @ before processing the first event sample
    void ProcessBegin();

    /// Function to evaluate input showers and determine a score
    bool Analyze(const EventData &data, const ParticleGraph &ps);

    /// Called after processing the last event sample
    void ProcessEnd(TFile* fout=nullptr);

  };
}
#endif

/** @} */ // end of doxygen group 
