#ifndef ERTOOL_ERANAJEREMYPI0_CXX
#define ERTOOL_ERANAJEREMYPI0_CXX

#include "ERAnaJeremyPi0.h"

namespace ertool {

  ERAnaJeremyPi0::ERAnaJeremyPi0(const std::string& name) : AnaBase(name)
  {}

  void ERAnaJeremyPi0::Reset()
  {}

  void ERAnaJeremyPi0::AcceptPSet(const ::fcllite::PSet& cfg)
  {}

  void ERAnaJeremyPi0::ProcessBegin()
  {
    _nTruePi0 = _nRecoPi0 = 0;
    
    // instantiate histograms for analysis
    _pi0True        = new TH1I("pi0True", "pi0 reco efficiency", 10, 0, 1);
    _pi0EnergyEff   = new TH1D("pi0EnergyEff", "pi0 energy reconstruction efficiency", 40, 0 ,2);
  }

  bool ERAnaJeremyPi0::Analyze(const EventData &data, const ParticleGraph &ps)
  {
    auto mc = MCParticleGraph();
    
    auto true_nodes = mc.GetParticleNodes(RecoType_t::kInvisible, 0, 111);
    auto reco_nodes = ps.GetParticleNodes(RecoType_t::kInvisible, 0, 111);
    
    _nTruePi0 += true_nodes.size();
    _nRecoPi0 += reco_nodes.size();
 
    if (true_nodes.size() == 1 && reco_nodes.size() == 1)
    {
      auto real_pi0 = ps.GetParticle(reco_nodes[0]);
      auto reco_pi0 = mc.GetParticle(true_nodes[0]);
      double energyEff = reco_pi0.Energy() / real_pi0.Energy();
      _pi0EnergyEff->Fill(energyEff);
    }
    return true;
  }

  void ERAnaJeremyPi0::ProcessEnd(TFile* fout)
  {
    TFile* f = new TFile("pi0Eff.root", "UPDATE");

    if (!(f->Get("EnergyEff")))
    {
      TH1* RecoEff = new TH1D("RecoEff", "% of pi0s reconstructed", 10, 0, 10);
      RecoEff->Write();
      delete RecoEff;
    }
    if (!(f->Get("EnergyEff")))
    {
      TH2* EnergyEff = new TH2D("EnergyEff", "pi0 energy reconstruction efficiency", 10, 0, 10, 40, 0, 2);
      EnergyEff->Write();
      delete EnergyEff;
    }
  
    
    TH1* h1 = (TH1*)f->Get("RecoEff");
    double eff = 100 * _nRecoPi0 / _nTruePi0;
    h1->SetBinContent(_bin, eff);
    h1->Write();
    delete h1;
    
    TH2* h2 = (TH2*)f->Get("EnergyEff");
    for (int i = 1; i < 41; i++)
    {
      double temp = _pi0EnergyEff->GetBinContent(i);
      h2->SetBinContent(_bin, i, temp);
    }
    h2->Write();
    delete h2;
    f->Close();
    
    delete _pi0True;
    delete _pi0EnergyEff;
    
    // draw all the histograms and save them to .png files
    //TCanvas *c = new TCanvas("c", "Canvas", 1000, 500);
    
    //_pi0True->Draw();
    //c->SaveAs("pi0eff.png");
    
    //_pi0EnergyEff->Draw();
    //c->SaveAs(Form("pi0energyeff_%i.png", _bin));
    
    //delete c;
  }
  

  // ******************** //
  // **** CHECK PI0S **** //
  // ******************** //
  /*
  void ERAlgoJeremyPi0::CheckPi0s(ParticleGraph graph)
  {
    
    std::cout << "[" << __FUNCTION__ << "] PI0 EVENT SUMMARY" << std::endl;
    std::cout << std::endl;
    
    auto pi0_ids = graph.GetParticleNodes(RecoType_t::kInvisible,0,111);
    for (auto const& pi0_id : pi0_ids)
    {
      auto pi0 = graph.GetParticle(pi0_id);
      std::cout << "[" << __FUNCTION__ << "] Node ID   " << pi0_id << std::endl;
      std::cout << "[" << __FUNCTION__ << "] Mass      " << pi0.Mass() << "MeV" << std::endl;
      std::cout << "[" << __FUNCTION__ << "] Energy    " << pi0.Energy() << "MeV" << std::endl;
      std::cout << "[" << __FUNCTION__ << "] Vertex    (" << pi0.Vertex()[0] << "," << pi0.Vertex()[1] << "," << pi0.Vertex()[2] << ")" << std::endl;
      std::cout << "[" << __FUNCTION__ << "] Momentum  (" << pi0.Momentum()[0] << "," << pi0.Momentum()[1] << "," << pi0.Momentum()[2] << ")" << std::endl;
      std::cout << std::endl;
    }
  }*/
  
}

#endif
