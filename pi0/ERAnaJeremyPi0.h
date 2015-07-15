/**
 * \file ERAnaJeremyPi0.h
 *
 * \ingroup Selection
 * 
 * \brief Class def header for a class ERAnaJeremyPi0
 *
 * @author jhewes15
 */

/** \addtogroup Selection

    @{*/

#ifndef ERTOOL_ERANAJEREMYPI0_H
#define ERTOOL_ERANAJEREMYPI0_H

#include "ERTool/Base/AnaBase.h"
#include "ERTool/Base/EventData.h"
#include "ERTool/Base/ParticleGraph.h"
#include <TH1.h>
#include <TCanvas.h>

namespace ertool {

  /**
     \class ERAnaJeremyPi0
     User custom Analysis class made by kazuhiro
   */
  class ERAnaJeremyPi0 : public AnaBase {
  
  public:

    /// Default constructor
    ERAnaJeremyPi0(const std::string& name="ERAnaJeremyPi0");

    /// Default destructor
    virtual ~ERAnaJeremyPi0(){}

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
    
  protected:
    int _nTruePi0;
    int _nRecoPi0;
    
    // histograms
    TH1* _pi0True;
    TH1* _pi0Reco;
    TH1* _pi0EnergyEff;
  };
}
#endif

/** @} */ // end of doxygen group 
