#ifndef ERTOOL_ERANANNBAR_CXX
#define ERTOOL_ERANANNBAR_CXX

#include "ERAnaNNbar.h"

namespace ertool {

  ERAnaNNbar::ERAnaNNbar(const std::string& name) : AnaBase(name)
  {}

  void ERAnaNNbar::Reset()
  {}

  void ERAnaNNbar::AcceptPSet(const ::fcllite::PSet& cfg)
  {}

  void ERAnaNNbar::ProcessBegin()
  {}

  bool ERAnaNNbar::Analyze(const EventData &data, const ParticleGraph &ps)
  { return true; }

  void ERAnaNNbar::ProcessEnd(TFile* fout)
  {}

}

#endif
