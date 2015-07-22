#ifndef NGAMMABASE_CXX
#define NGAMMABASE_CXX

#include "NGammaBase.h"

namespace ertool {
  
  
  // ********************* //
  // **** CONSTRUCTOR **** //
  // ********************* //

  NGammaBase::NGammaBase()
    : _verbose(false)
    , _pset("NGamma")
  {
    // set the external variables to zero
    _extEnergy = 0.;
    _extMomentum = new geoalgo::Vector(0.,0.,0.);
    
    // instantiate all the roofit vars, pdfs & datasets
    _energyVar     = new RooRealVar("energyVar", "Energy [MeV] variable", 0.1, 2000.);
    _energyPdf     = _factory.Gaus("energyPdf", *_energyVar);
    _energyData    = new RooDataSet("energyData", "NGamma energy data for PDF training",RooArgSet(*_energyVar));
    
    _momentumVar   = new RooRealVar("momentumVar", "Momentum [MeV] variable", 0.1, 2000.);
    _momentumPdf   = _factory.Gaus("momentumPdf", *_momentumVar);
    _momentumData  = new RooDataSet("momentumData", "NGamma momentum data for PDF training", RooArgSet(*_momentumVar));
    
    _massVar       = new RooRealVar("massVar", "Invariant mass [MeV] variable", 0.1, 2000.);
    _massPdf       = _factory.Gaus("massPdf", *_massVar);
    _massData      = new RooDataSet("massData", "NGamma invariant mass data for PDF training", RooArgSet(*_massVar));
  }
  
  
  // ********************** //
  // **** ACCEPT P SET **** //
  // ********************** //
  
  void NGammaBase::AcceptPSet(const ::fcllite::PSet& cfg)
  {
    // temporary mean & variable holders
    //RooRealVar *meanVar, *sigmaVar;
    
    // pull out info for pdf instantiation
    auto p = cfg.get_pset("NGamma");
    
    // ENERGY PDF PARAMETERS
    if (p.contains_value("energy_params"))
    {
      _useEnergyPdf = true;
      auto darray = p.get<std::vector<double> >("energy_params");
      
      double mean  = darray[0];
      double sigma = darray[1];
      
      double min = mean - (3 * sigma);
      double max = mean + (3 * sigma);
      
      // set range of mass pdf gaussian parameter
      _energyVar->setRange(min, max);
      
      RooRealVar *meanVar, *sigmaVar;
      
      // set value of mass pdf gaussian mean
      meanVar = (RooRealVar*)(_energyPdf->getVariables()->find("energyPdf_Gaus_mean"));
      SetParam(meanVar, mean);
      
      // set value of energy pdf gaussian sigma
      sigmaVar = (RooRealVar*)(_energyPdf->getVariables()->find("energyPdf_Gaus_sigma"));
      SetParam(sigmaVar, sigma);
      
      if (_verbose)
        std::cout << "[" << __FUNCTION__ << "] Set new parameters for energy PDF: mean " << meanVar->getVal() << ", sigma " << sigmaVar->getVal() << ", range " << _energyVar->getMin() << " ==> " << _energyVar->getMax() << std::endl;
    }
    
    // MOMENTUM PDF PARAMETERS
    if (p.contains_value("momentum_params"))
    {
      _useMomentumPdf = true;
      auto darray = p.get<std::vector<double> >("momentum_params");
      
      // pull variables out of config file
      double mean  = darray[0];
      double sigma = darray[1];
      double min   = mean - (3 * sigma);
      double max   = mean + (3 * sigma);
      
      // input all the values into roofit
      _momentumVar->setRange(min, max);
      RooRealVar *meanVar, *sigmaVar;
      meanVar = (RooRealVar*)(_momentumPdf->getVariables()->find("momentumPdf_Gaus_mean"));
      SetParam(meanVar, mean);
      sigmaVar = (RooRealVar*)(_momentumPdf->getVariables()->find("momentumPdf_Gaus_sigma"));
      SetParam(sigmaVar, sigma);
      
      if (_verbose)
        std::cout << "[" << __FUNCTION__ << "] Set new parameters for momentum PDF: mean " << meanVar->getVal() << ", sigma " << sigmaVar->getVal() << ", range " << _momentumVar->getMin() << " ==> " << _momentumVar->getMax() << std::endl;
    }
    
    // MASS PDF PARAMETERS
    if (p.contains_value("mass_params"))
    {
      _useMassPdf = true;
      auto darray = p.get<std::vector<double> >("mass_params");
      
      double mean  = darray[0];
      double sigma = darray[1];
      double min   = mean - (3 * sigma);
      double max   = mean + (3 * sigma);
      
      _massVar->setRange(min, max);
      RooRealVar *meanVar, *sigmaVar;
      meanVar = (RooRealVar*)(_massPdf->getVariables()->find("massPdf_Gaus_mean"));
      SetParam(meanVar, mean);
      sigmaVar = (RooRealVar*)(_massPdf->getVariables()->find("massPdf_Gaus_sigma"));
      SetParam(sigmaVar, sigma);
      
      if (_verbose)
        std::cout << "[" << __FUNCTION__ << "] Set new parameters for mass PDF: mean " << meanVar->getVal() << ", sigma " << sigmaVar->getVal() << ", range " << _massVar->getMin() << " ==> " << _massVar->getMax() << std::endl;
    }
  }
  
  
  // ************************** //
  // **** GAMMA COMPARISON **** //
  // ************************** //
  
  CombinationScoreSet_t NGammaBase::GammaComparison(const ParticleGraph& graph, int n)
  {
    CombinationScoreSet_t scores;
    
    // if you don't specify what size combination you want...

    int a, b;
    if (n == 0)
    {
      a = 1;
      b = graph.GetParticleNodes(RecoType_t::kShower, 0, 22).size();
    }
    else if (n > graph.GetParticleNodes(RecoType_t::kShower, 0, 22).size())
      return scores;
    else
      a = b = n;
    
    // ...then the code automatically loops over all possible sizes of combination
    for (int i = a; i < b + 1; i++)
    {
      // get all possible combinations for a given combination size
      CombinationSet_t combinations = graph.GetNodeCombinations(i, RecoType_t::kShower, 0, 22);
        
      // loop over all possible combinations
      for (auto const& combination : combinations)
      {
          
        // perform a likelihood calculation for that combination
        double l = Likelihood(graph, combination);
        
        // add up all the energy in that combination (& spit out some output if desired)
        if (_verbose)
        {
          std::cout << "[" << __FUNCTION__ << "] Comparing showers";
          for (auto const& particle : combination)
            std::cout << " " << graph.GetParticle(particle).ID();
          std::cout << " with score " << l << std::endl;
        }
        
        // add the combination & associated score to the "combination score set" object
        CombinationScore_t score(combination, l);
        scores.push_back(score);
      }
    }
    return scores;
  }

  
  // ************************* //
  // **** EVENT SELECTION **** //
  // ************************* //
  
  CombinationScoreSet_t NGammaBase::EventSelection(CombinationScoreSet_t candidates, double threshold)
  {
    // sort combinations in descending score order & remove any below a certain threshold
    std::sort(candidates.begin(), candidates.end(), [](CombinationScore_t a, CombinationScore_t b){ return (a.second > b.second); });
    candidates.erase(std::remove_if(candidates.begin(), candidates.end(), [=](CombinationScore_t a){ return (a.second < threshold); }), candidates.end());
    
    // create a new event object
    CombinationScoreSet_t event;
    
    // loop over all potential candidates
    while (!candidates.empty())
    {
      // grab the candidate with the best score & add to event
      Combination_t gammas = candidates.front().first;
      event.push_back(candidates.front());
      
      // instantiate a temporary object from your candidates
      CombinationScoreSet_t temp = candidates;
      candidates.clear();
      
      // loop back over your temporary object, saving only candidates that don't clash with the just-selected combination
      for (auto const& candidate : temp)
      {
        bool overlap = false;
        for (auto const& id1 : candidate.first)
          for (auto const& id2 : gammas)
            if (id1 == id2)
              overlap = true;
        if (!overlap)
          candidates.push_back(candidate);
      }
      
      // clear your temporary object
      temp.clear();
    }
    return event;
  }

  
  // i've been staring at code for too long
  // here's a mixtape
  
  // winter of the electric beach -- olan mill
  // reading in bed -- emily haines & the soft skeleton
  // loma linda -- nik freitas
  // i -- kendrick lamar
  // candy's room -- lifter puller
  // nervous -- wavves x cloud nothings
  // happy to see me -- hop along
  // psalm -- m. ward
  // sunrise at wildflower hill retirement center -- kaki king
  // car 24 -- seatbelts
  // transcendental youth -- the mountain goats
  
  // ok, back to coding now
  
  
  // ******************** //
  // **** LIKELIHOOD **** //
  // ******************** //
  
  double NGammaBase::Likelihood(ParticleGraph graph, Combination_t combination)
  {
    // first check to make sure at least one pdf is instantiated
    if(_useEnergyPdf == false && _useMomentumPdf == false && _useMassPdf == false)
    {
      std::cout << "Error! Trying to perform a likelihood calculation, but no PDFs have been instantiated!" << std::endl;
      return 0;
    }
    
    double l = 1;
    
    // run through all the pdfs, and contribute to the score if the pdf is used

    if (_useEnergyPdf)
    {
      // calculate energy of combination
      double energy = _extEnergy;
      for (auto const& particle : combination)
        energy += graph.GetParticle(particle).Energy();
      
      // get likelihood from pdf
      _energyVar->setVal(energy);
      double e = _energyPdf->getVal(*_energyVar);
      double emax = 1 / _energyPdf->getNorm(*_energyVar);
      double p = e / emax;
      l *= p;
    }
    
    if (_useMomentumPdf)
    {
      // calculate net momentum of combination
      geoalgo::Vector momentum = *_extMomentum;
      for (auto const& particle : combination)
        momentum += graph.GetParticle(particle).Momentum();
      
      // get likelihood from pdf
      _momentumVar->setVal(momentum.Length());
      double e = _momentumPdf->getVal(*_momentumVar);
      double emax = 1 / _momentumPdf->getNorm(*_momentumVar);
      double p = e / emax;
      l *= p;
    }
    
    if (_useMassPdf)
    {
      // make sure we're dealing with a pair
      if (combination.size() != 2)
      {
        std::cout << "[" << __FUNCTION__ << "] Error using mass PDF for likelihood calculation!" << std::endl;
        std::cout << "[" << __FUNCTION__ << "] Mass PDF is designed specifically for calculating the invariant mass of a pair of showers." << std::endl;
        std::cout << "[" << __FUNCTION__ << "] The code is trying to run for a combination of " << combination.size() << " showers." << std::endl;
        return 0;
      }
      
      // calculate invariant mass of pair
      auto particle1 = graph.GetParticle(combination[0]);
      auto particle2 = graph.GetParticle(combination[1]);
      double angle   = particle1.Momentum().Angle(particle2.Momentum());
      double invMass = sqrt(2 * particle1.Energy() * particle2.Energy() * (1 - cos(angle)));
      
      // get likelihood from pdf
      _massVar->setVal(invMass);
      double e = _massPdf->getVal(*_massVar);
      double emax = 1 / _massPdf->getNorm(*_massVar);
      double p = e / emax;
      l *= p;
    }
    return l;
  }
  
  // https://twitter.com/CuteEmergency/status/619165430582804480
  
  
  // *************************************** //
  // **** ROOFIT DATA FILLING FUNCTIONS **** //
  // *************************************** //
  
  void NGammaBase::FillEnergyPdf(double energy)
  {
    // make sure energy value is non-negative
    if (energy >= 0)
    {
      // tell the code that you're using this pdf
      _useEnergyPdf = true;
      
      // set roofit variable & add to data set
      _energyVar->setVal(energy);
      _energyData->add(RooArgSet(*_energyVar));
    }
    
    // complain if the energy is negative
    else
      if (_verbose)
        std::cout << "[" << __FUNCTION__ << "] Energy value (" << energy << ") cannot be negative!" << std::endl;
  }
  
  void NGammaBase::FillMomentumPdf(geoalgo::Vector momentum)
  {
    // tell the code you're using this pdf
    _useMomentumPdf = true;
    
    // set roofit variable & add to data set
    _momentumVar->setVal(momentum.Length());
    _momentumData->add(RooArgSet(*_momentumVar));
  }
  
  void NGammaBase::FillMassPdf(double mass)
  {
     // tell the code you're using this pdf
    _useMassPdf = true;
    
    // set roofit variable & add to data set
    _massVar->setVal(mass);
    _massData->add(RooArgSet(*_massVar));
  }
  
  fcllite::PSet NGammaBase::GetParams()
  {
    // make sure at least one pdf has been instantiated before opening output file
    if (!(_useEnergyPdf || _useMomentumPdf || _useMassPdf) && _verbose)
      std::cout << "[" << __FUNCTION__ << "] No PDFs were instantiated! Not saving output." << std::endl;
    
    // pointers for mean & sigma variables
    RooRealVar *meanVar, *sigmaVar;
    
    // save energy params
    if (_useEnergyPdf)
    {
      if (_verbose)
      {
        std::cout << "[" << __FUNCTION__ << "] ";
        _energyData->Print();
      }
      
      // get pointers to mean & sigma
      meanVar = (RooRealVar*)(_energyPdf->getVariables()->find("energyPdf_Gaus_mean"));
      sigmaVar = (RooRealVar*)(_energyPdf->getVariables()->find("energyPdf_Gaus_sigma"));
      
      // stop it from crashing somehow
      double min, max, mean;
      _energyData->getRange(*_energyVar, min, max);
      min *= 0.9;
      max *= 1.1;
      mean = _energyData->mean(*_energyVar);
      meanVar->setRange(5, max);
      sigmaVar->setRange(5, max);
      meanVar->setVal(mean);
      sigmaVar->setVal(mean);
      
      // fit pdf to data
      RooFitResult* energyResult = _energyPdf->fitTo(*_energyData, RooFit::Save(), RooFit::PrintLevel(-1));
      if (_verbose)
        energyResult->Print();
      
      // plot pdf & data
      TCanvas *c = new TCanvas("c","",1000,500);
      RooPlot* frame = _energyVar->frame();
      _energyData->plotOn(frame, RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack));
      _energyPdf->plotOn(frame, RooFit::LineColor(kRed));
      frame->Draw();
      c->SaveAs("./energyPdf.png");
      delete c;

      // save the parameter set to a config file
      std::vector<double> darray;
      auto& params = OutputPSet();
      
      darray.push_back(meanVar->getVal());
      darray.push_back(sigmaVar->getVal());
      
      params.add_value("energy_params",::fcllite::VecToString(darray));
    }
    
    // save momentum params
    if (_useMomentumPdf)
    {
      if (_verbose)
      {
        std::cout << "[" << __FUNCTION__ << "] ";
        _momentumData->Print();
      }
      
      // get pointers to mean & sigma
      meanVar = (RooRealVar*)(_momentumPdf->getVariables()->find("momentumPdf_Gaus_mean"));
      sigmaVar = (RooRealVar*)(_momentumPdf->getVariables()->find("momentumPdf_Gaus_sigma"));
      
      // do, uh, something
      double min, max, mean;
      _momentumData->getRange(*_momentumVar, min, max);
      std::cout << "Minimum is " << min << ", maximum is " << max;
      min *= 0.9;
      max *= 1.1;
      mean = _momentumData->mean(*_momentumVar);
      meanVar->setRange(5, max);
      sigmaVar->setRange(5, max);
      meanVar->setVal(mean);
      sigmaVar->setVal(mean);
      
      // fit pdf to data
      RooFitResult* momentumResult = _momentumPdf->fitTo(*_momentumData, RooFit::Save(), RooFit::PrintLevel(-1));
      if (_verbose)
        momentumResult->Print();
      
      // plot pdf & data
      TCanvas *c = new TCanvas("c","",1000,500);
      RooPlot* frame = _momentumVar->frame();
      _momentumData->plotOn(frame, RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack));
      _momentumPdf->plotOn(frame, RooFit::LineColor(kRed));
      frame->Draw();
      c->SaveAs("momentumPdf.png");
      delete c;
      
      // save the parameter set to a config file
      std::vector<double> darray;
      auto& params = OutputPSet();
      
      darray.push_back(meanVar->getVal());
      darray.push_back(sigmaVar->getVal());
      
      params.add_value("momentum_params",::fcllite::VecToString(darray));
    }
    
    // save mass params
    if (_useMassPdf)
    {
      if (_verbose)
      {
        std::cout << "[" << __FUNCTION__ << "] ";
        _massData->Print();
      }
      
      // get pointers to mean & sigma
      meanVar = (RooRealVar*)(_massPdf->getVariables()->find("massPdf_Gaus_mean"));
      sigmaVar = (RooRealVar*)(_massPdf->getVariables()->find("massPdf_Gaus_sigma"));
      
      // set ranges of variables (kind of a black art tbh)
      double min, max, mean;
      _massData->getRange(*_massVar, min, max);
      min *= 0.9;
      max *= 1.1;
      mean = _massData->mean(*_massVar);
      _massVar->setRange(0.1, max);
      meanVar->setRange(5, max);
      sigmaVar->setRange(5, max);
      meanVar->setVal(mean);
      sigmaVar->setVal(mean);
      
      // fit pdf to data
      RooFitResult* massResult = _massPdf->fitTo(*_massData, RooFit::Save(), RooFit::PrintLevel(-1));
      if (_verbose)
        massResult->Print();
      
      // plot pdf & data
      TCanvas *c = new TCanvas("c","",1000,500);
      RooPlot* frame = _massVar->frame();
      _massData->plotOn(frame, RooFit::MarkerColor(kBlack), RooFit::LineColor(kBlack));
      _massPdf->plotOn(frame, RooFit::LineColor(kRed));
      frame->GetXaxis()->SetLimits(0,250);
      frame->Draw();
      c->SaveAs("./massPdf.png");
      delete c;
      
      // save the parameter set to a config file
      std::vector<double> darray;
      auto& params = OutputPSet();
      
      darray.push_back(meanVar->getVal());
      darray.push_back(sigmaVar->getVal());
      
      params.add_value("mass_params",::fcllite::VecToString(darray));
    }
    return OutputPSet();
  }
  
  
  // ******************* //
  // **** SET PARAM **** //
  // ******************* //
  
  void NGammaBase::SetParam(RooRealVar* var, double val)
  {
    double max = var->getMax();
    double min = var->getMin();
    
    if ( val > max )
      var->setMax(1.1 * val);
    else if (val < min)
      var->setMin(0.9 * val);
    
    var->setVal(val);
  }
  
}

#endif
