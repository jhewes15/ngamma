#ifndef ERTOOL_ERALGOJEREMYPI0_CXX
#define ERTOOL_ERALGOJEREMYPI0_CXX

#include "ERAlgoJeremyPi0.h"

namespace ertool {

  ERAlgoJeremyPi0::ERAlgoJeremyPi0(const std::string& name) : AlgoBase(name)
  {
    _nGamma = new NGammaBase("pi0");
  }

  void ERAlgoJeremyPi0::Reset()
  {}

  void ERAlgoJeremyPi0::AcceptPSet(const ::fcllite::PSet& cfg)
  {
    if (!_trainingMode)
      _nGamma->AcceptPSet(cfg);
  }

  void ERAlgoJeremyPi0::ProcessBegin()
  {}

  bool ERAlgoJeremyPi0::Reconstruct(const EventData &data, ParticleGraph& graph)
  {
    if (_trainingMode)
    {
      auto nodes = graph.GetParticleNodes(RecoType_t::kShower,0,22);
      if (nodes.size() == 2)
      {
        auto particle1 = graph.GetParticle(nodes[0]);
        auto particle2 = graph.GetParticle(nodes[1]);
        double angle = particle1.Momentum().Angle(particle2.Momentum());
        double invMass = sqrt(2 * particle1.Energy() * particle2.Energy() * (1 - cos(angle)));
        _nGamma->FillMassPdf(invMass);
      }
    }
    
    else
    {
      CombinationScoreSet_t candidates = _nGamma->GammaComparison(graph, 2);
      CombinationScoreSet_t event = _nGamma->EventSelection(candidates, 0.7);
      AddPi0s(graph, event);
      CheckPi0s(graph);
    }
    
    return true;
  }

  void ERAlgoJeremyPi0::ProcessEnd(TFile* fout)
  {
    if (_trainingMode)
      _nGamma->SaveParams();
  }
  
  // ****************** //
  // **** ADD PI0S **** //
  // ****************** //
  
  void ERAlgoJeremyPi0::AddPi0s(ParticleGraph &graph, CombinationScoreSet_t particles)
  {
    // this comes in useful later:
    int pdg_pi0 = 111;
    
    // loop over each identified pi0
    for (auto const& particle : particles)
    {
      // get ids of two gamma showers
      double gamma1_id = particle.first[0];
      double gamma2_id = particle.first[1];
      
      // get pointers to two gamma showers
      Particle *gamma1 = &(graph.GetParticle(gamma1_id));
      Particle *gamma2 = &(graph.GetParticle(gamma2_id));
      
      // identify them as siblings (which creates "parent" object if it doesn't exist already)
      graph.SetSiblings(gamma1_id, gamma2_id, particle.second);
      
      // get id of parent object
      double parent_id = gamma1->Parent();
      
      // get a pointer to the parent object
      Particle *parent = &(graph.GetParticle(parent_id));
      
      // figure out the pi0's information
      geoalgo::Vector X (( gamma1->Vertex()[0] + gamma2->Vertex()[0]) /2, (gamma1->Vertex()[1] + gamma2->Vertex()[1]) /2, (gamma1->Vertex()[2] + gamma2->Vertex()[2]) /2 );
      geoalgo::Vector P ( gamma1->Momentum()[0] + gamma2->Momentum()[0], gamma1->Momentum()[1] + gamma2->Momentum()[1], gamma1->Momentum()[2] + gamma2->Momentum()[2] );
      parent->SetParticleInfo(pdg_pi0, ParticleMass(pdg_pi0), X, P, particle.second);
    }
  }
  
  
  // ******************** //
  // **** CHECK PI0S **** //
  // ******************** //
  
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
  }
  
  // ******************************** //
  // **** GAMMA SHOWER GENERATOR **** //
  // ******************************** //
  
  Shower ERAlgoJeremyPi0::Generate()
  {
    // instantiate random number generator
    TRandom3 rand;
    
    rand.SetSeed(0);
    
    // randomly generate position & momentum
    double x = rand.Gaus(0,5);
    double y = rand.Gaus(0,5);
    double z = rand.Gaus(0,5);
    
    double px = rand.Gaus(0,5);
    double py = rand.Gaus(0,5);
    double pz = rand.Gaus(0,5);
    
    // if possible, get de/dx using empart params
    
    double dedx, energy;
    
    if (_empartMean1 == -1) dedx = 5;
    else
    {
      if (rand.Integer(2) == 0)
        dedx = rand.Gaus(_empartMean1,_empartSigma1);
      else
        dedx = rand.Gaus(_empartMean2,_empartSigma2);
    }
    
    
    // generate some energy yo
    /*
     if (_useEnergyPdf)
     {
     RooRealVar *meanVar, *sigmaVar;
     meanVar  = (RooRealVar*)(_energyPdf->getVariables()->find("energyPdf_Gaus_mean"));
     double mean  = meanVar->getVal();
     sigmaVar = (RooRealVar*)(_energyPdf->getVariables()->find("energyPdf_Gaus_sigma"));
     double sigma = sigmaVar->getVal();
     energy = rand.Gaus(mean,sigma);
     }
     else
     energy = 100;
     */
    
    energy = rand.Uniform(0,1500);
    
    // instantiate the shower object
    Shower s( geoalgo::Vector(x,y,z), geoalgo::Vector(px,py,pz), 30, 15);
    s._energy = energy;
    s._dedx = dedx;
    
    if (_verbose)
      std::cout << "[" << __FUNCTION__ << "] Gamma generated with position (" << x << "," << y << "," << z << "), momentum (" << px << "," << py << "," << pz << "), dedx " << dedx << ", energy " << energy << std::endl;
    
    return s;
  }
  
}

#endif
