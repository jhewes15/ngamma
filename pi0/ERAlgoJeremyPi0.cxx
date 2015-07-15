#ifndef ERTOOL_ERALGOJEREMYPI0_CXX
#define ERTOOL_ERALGOJEREMYPI0_CXX

#include "ERAlgoJeremyPi0.h"

namespace ertool {

  ERAlgoJeremyPi0::ERAlgoJeremyPi0(const std::string& name) : AlgoBase(name)
  {
    _name = name;
    _nGamma = new NGammaBase();
  }

  void ERAlgoJeremyPi0::Reset()
  {}

  void ERAlgoJeremyPi0::AcceptPSet(const ::fcllite::PSet& cfg)
  {
    // get external parameters if not training PDFs
    if (!_trainingMode)
    {
      auto ngamma_pset = cfg.get_pset(_name);
      _nGamma->AcceptPSet(ngamma_pset);
    }
  }

  void ERAlgoJeremyPi0::ProcessBegin()
  {}

  bool ERAlgoJeremyPi0::Reconstruct(const EventData &data, ParticleGraph& graph)
  {
    // first case: training pdfs!
    if (_trainingMode)
    {
      auto nodes = graph.GetParticleNodes(RecoType_t::kShower,0,22);
      if (nodes.size() == 2)
      {
        // get the two gamma shower particles
        auto particle1 = graph.GetParticle(nodes[0]);
        auto particle2 = graph.GetParticle(nodes[1]);
        
        // calculate angle between showers
        double angle = particle1.Momentum().Angle(particle2.Momentum());
        
        // calculate invariant mass
        double invMass = sqrt(2 * particle1.Energy() * particle2.Energy() * (1 - cos(angle)));
        
        // make sure energies are above a certain threshold
        if (invMass < 300)
        {
          // add to PDF
          _nGamma->FillMassPdf(invMass);
        }
      }
    }
    
    // second case: using pdfs to make selection!
    else
    {
      // score all possible combinations of gammas
      CombinationScoreSet_t candidates = _nGamma->GammaComparison(graph, 2);
      
      // select the best candidates as pi0s
      CombinationScoreSet_t event = _nGamma->EventSelection(candidates, 0.7);
      
      // add the pi0s to the particle graph, then check it worked properly
      AddPi0s(graph, event);
    }
    
    return true;
  }

  void ERAlgoJeremyPi0::ProcessEnd(TFile* fout)
  {
    // if in training mode, save the trained parameters for future use
    if (_trainingMode)
    {
      auto& params = OutputPSet();
      params.add_pset(_nGamma->GetParams());
    }
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
