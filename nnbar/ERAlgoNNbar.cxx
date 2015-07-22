#ifndef ERTOOL_ERALGONNBAR_CXX
#define ERTOOL_ERALGONNBAR_CXX

#include "ERAlgoNNbar.h"

namespace ertool {

  ERAlgoNNbar::ERAlgoNNbar(const std::string& name) : AlgoBase(name)
  {
    _name = name;
    _nGamma = new NGammaBase();
  }

  void ERAlgoNNbar::Reset()
  {}

  void ERAlgoNNbar::AcceptPSet(const ::fcllite::PSet& cfg)
  {
    // get external pdf params if in selection mode
    if (!_trainingMode)
    {
      auto ngamma_pset = cfg.get_pset(_name);
      _nGamma->AcceptPSet(ngamma_pset);
    }
  }

  void ERAlgoNNbar::ProcessBegin()
  {}

  bool ERAlgoNNbar::Reconstruct(const EventData &data, ParticleGraph& graph)
  {
    // for pdf training
    if (_trainingMode)
    {
      // instantiate variables
      double energy = 0;
      geoalgo::Vector momentum(0,0,0);
      
      // get final state particles
      auto gammaNodes = graph.GetParticleNodes(RecoType_t::kShower, 0, 22);
      auto piPNodes   = graph.GetParticleNodes(RecoType_t::kTrack, 0, 211);
      auto piMNodes   = graph.GetParticleNodes(RecoType_t::kTrack, 0, -211);
      
      // add up energy & momentum from all of them
      for (auto const& gammaNode : gammaNodes)
      {
        auto gamma = graph.GetParticle(gammaNode);
        energy    += gamma.Energy();
        momentum  += gamma.Momentum();
      }
      for (auto const& piPNode : piPNodes)
      {
        auto piP   = graph.GetParticle(piPNode);
        energy    += piP.Energy();
        momentum  += piP.Momentum();
      }
      for (auto const& piMNode : piMNodes)
      {
        auto piM   = graph.GetParticle(piMNode);
        energy    += piM.Energy();
        momentum  += piM.Momentum();
      }
      
      // add these to pdfs
      _nGamma->FillEnergyPdf(energy);
      _nGamma->FillMomentumPdf(momentum);
    }
    
    // for event selection
    else
    {
      // get number of pi plus & pi minus in event
      int nPiP = graph.GetParticleNodes(RecoType_t::kTrack, 0, 211).size();
      int nPiM = graph.GetParticleNodes(RecoType_t::kTrack, 0, -211).size();
      
      // loop over each combination of them
      for (int i = 0; i < nPiP; i++)
      {
        for (int j = 0; j < nPiM; j++)
        {
          auto piPCombinations = graph.GetNodeCombinations(i, RecoType_t::kTrack, 0, 211);
          auto piMCombinations = graph.GetNodeCombinations(j, RecoType_t::kTrack, 0, -211);
          for (auto const& piPCombination : piPCombinations)
          {
            for (auto const& piMCombination : piMCombinations)
            {
              // for each combination, add up total energy & momentum
              double energy = 0;
              geoalgo::Vector momentum(0,0,0);
              for (auto const& piPNode : piPCombination)
              {
                auto piP  = graph.GetParticle(piPNode);
                energy   += piP.Energy();
                momentum += piP.Momentum();
              }
              for (auto const& piMNode : piMCombination)
              {
                auto piM  = graph.GetParticle(piMNode);
                energy   += piM.Energy();
                momentum += piM.Momentum();
              }
              _nGamma->SetExternalEnergy(energy);
              _nGamma->SetExternalMomentum(momentum);
              
              auto candidates = _nGamma->GammaComparison(graph);
              auto events     = _nGamma->EventSelection(candidates,0.9);
              
              if (events.size() > 0)
                std::cout << "[" << __FUNCTION__ << "] Candidate found!" << std::endl;
            }
          }
        }
      }
    }
    return true;
  }

  void ERAlgoNNbar::ProcessEnd(TFile* fout)
  {
    // if in training mode, save trained parameters for future use
    if (_trainingMode)
    {
      auto& params = OutputPSet();
      params.add_pset(_nGamma->GetParams());
    }
  }

}

#endif
