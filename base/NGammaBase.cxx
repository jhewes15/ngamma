#ifndef NGAMMABASE_CXX
#define NGAMMABASE_CXX

#include "NGammaBase.h"

namespace ertool {
  
  
  // ********************* //
  // **** CONSTRUCTOR **** //
  // ********************* //

  NGammaBase::NGammaBase(std::string name)
    : _verbose(false)
    , _pset(name)
  {
    // set the external variables to zero
    _extEnergy = 0.;
    _extMomentum = new geoalgo::Vector(0.,0.,0.);
    
    // instantiate all the roofit vars, pdfs & datasets
    _energyVar     = new RooRealVar("energyVar", "Energy [MeV] variable", 0.1, 10000.);
    _energyPdf     = _factory.Gaus("energyPdf", *_energyVar);
    _energyData    = new RooDataSet("energyData", "NGamma energy data for PDF training",RooArgSet(*_energyVar));
    
    _momentumVar   = new RooRealVar("momentumVar", "Momentum [MeV] variable", 0.1, 10000.);
    _momentumPdf   = _factory.Gaus("momentumPdf", *_momentumVar);
    _momentumData  = new RooDataSet("momentumData", "NGamma momentum data for PDF training", RooArgSet(*_momentumVar));
    
    _massVar       = new RooRealVar("massVar", "Invariant mass [MeV] variable", 0.1, 10000.);
    _massPdf       = _factory.Gaus("massPdf", *_massVar);
    _massData      = new RooDataSet("massData", "NGamma invariant mass data for PDF training", RooArgSet(*_massVar));
    
    _radLenVar     = new RooRealVar("radLenVar", "Radiation length [cm] variable", 0.1, 100.);
    _radLenPdf     = _factory.Gaus("radLenPdf", *_radLenVar);
    _radLenData    = new RooDataSet("radLenData", "NGamma radiation length data for PDF training", RooArgSet(*_radLenVar));
  }
  
  
  // ********************** //
  // **** ACCEPT P SET **** //
  // ********************** //
  
  void NGammaBase::AcceptPSet(const ::fcllite::PSet& cfg)
  {
    // temporary mean & variable holders
    RooRealVar *meanVar, *sigmaVar;
    
    // pull out info for pdf instantiation
    auto p = cfg.get_pset("NGamma");
    
    // ENERGY PDF PARAMETERS
    if (p.contains_value("mass_params"))
    {
      _useEnergyPdf = true;
      auto darray = p.get<std::vector<double> >("mass_params");
      
      double mean  = darray[0];
      double sigma = darray[1];
      
      double min = mean - (5 * sigma);
      double max = mean + (5 * sigma);
      
      // set range of mass pdf gaussian parameter
      _energyVar->setRange(min, max);
      
      // set value of mass pdf gaussian mean
      meanVar = (RooRealVar*)(_massPdf->getVariables()->find("massPdf_Gaus_mean"));
      //meanVar->setRange(0.9 * mean, 1.1 * mean);
      meanVar->setVal(mean);
      
      // set value of energy pdf gaussian sigma
      sigmaVar = (RooRealVar*)(_massPdf->getVariables()->find("massPdf_Gaus_sigma"));
      //sigmaVar->setRange(0.9 * sigma, 1.1 * sigma);
      sigmaVar->setVal(sigma);
      
      if (_verbose)
        std::cout << "[" << __FUNCTION__ << "] Set new parameters for mass PDF: mean " << meanVar->getVal() << ", sigma " << sigmaVar->getVal() << ", range " << _massVar->getMin() << " ==> " << _massVar->getMax() << std::endl;
    }
    
    // MOMENTUM PDF PARAMETERS
    if (p.contains_value("momentum_params"))
    {
      _useMomentumPdf = true;
      auto darray = p.get<std::vector<double> >("momentum_params");
      
      // pull variables out of config file
      double mean  = darray[0];
      double sigma = darray[1];
      double min   = mean - (5 * sigma);
      double max   = mean + (5 * sigma);
      
      // input all the values into roofit
      _momentumVar->setRange(min, max);
      meanVar = (RooRealVar*)(_momentumPdf->getVariables()->find("momentumPdf_Gaus_mean"));
      //meanVar->setRange(0.9 * xmean, 1.1 * xmean);
      meanVar->setVal(mean);
      sigmaVar = (RooRealVar*)(_momentumPdf->getVariables()->find("momentumPdf_Gaus_sigma"));
      //sigmaVar->setRange(0.9 * xsigma, 1.1 * xsigma);
      sigmaVar->setVal(sigma);
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
      double min   = mean - (5 * sigma);
      double max   = mean + (5 * sigma);
      
      _massVar->setRange(min, max);
      
      meanVar = (RooRealVar*)(_massPdf->getVariables()->find("massPdf_Gaus_mean"));
      meanVar->setVal(mean);
      
      sigmaVar = (RooRealVar*)(_massPdf->getVariables()->find("massPdf_Gaus_sigma"));
      sigmaVar->setVal(sigma);
      
      if (_verbose)
        std::cout << "[" << __FUNCTION__ << "] Set new parameters for mass PDF: mean " << meanVar->getVal() << ", sigma " << sigmaVar->getVal() << ", range " << _massVar->getMin() << " ==> " << _massVar->getMax() << std::endl;
    }
  
    // RADIATION LENGTH PDF PARAMETERS
    if (p.contains_value("radlen_params"))
    {
      _useRadLenPdf = true;
      auto darray = p.get<std::vector<double> >("radlen_params");
    
      double mean  = darray[0];
      double sigma = darray[1];
      double min   = mean - (5 * sigma);
      double max   = mean + (5 * sigma);
    
      _radLenVar->setRange(min, max);
    
      meanVar = (RooRealVar*)(_radLenPdf->getVariables()->find("radLenPdf_Gaus_mean"));
      meanVar->setVal(mean);
    
      sigmaVar = (RooRealVar*)(_radLenPdf->getVariables()->find("radLenPdf_Gaus_sigma"));
      sigmaVar->setVal(sigma);
    
      if (_verbose)
        std::cout << "[" << __FUNCTION__ << "] Set new parameters for radiation length PDF: mean " << meanVar->getVal() << ", sigma " << sigmaVar->getVal() << ", range " << _radLenVar->getMin() << " ==> " << _radLenVar->getMax() << std::endl;
    }
  }
  
  
  // ************************** //
  // **** GAMMA COMPARISON **** //
  // ************************** //
  
  CombinationScoreSet_t NGammaBase::GammaComparison(const ParticleGraph& graph, int n)
  {
    CombinationScoreSet_t scores;
    
    // if you don't specify what size combination you want...
    if (n == 0)
    {
      // ...then the code automatically loops over all possible sizes of combination
      for (int i = 0; i < graph.GetParticleNodes(RecoType_t::kShower,0,22).size(); i++)
      {
        
        // get all possible combinations for a given combination size
        CombinationSet_t combinations = graph.GetNodeCombinations(i+1, RecoType_t::kShower, 0, 22);
        
        // loop over all possible combinations
        for (auto const& combination : combinations)
        {
          // add up all the energy in that combination (& spit out some output if desired)
          if (_verbose)
          {
            std::cout << "[" << __FUNCTION__ << "] Comparing showers";
            for (auto const& particle : combination)
              std::cout << " " << graph.GetParticle(particle).ID();
          }
          
          // perform a likelihood calculation for that combination
          double l = Likelihood(graph, combination);
          
          // add the combination & associated score to the "combination score set" object
          CombinationScore_t score(combination, l);
          scores.push_back(score);
        }
      }
    }
    
    // does the same as above, but for one specific combination size
    else if (n > 0 && n <= graph.GetParticleNodes(RecoType_t::kShower,0,22).size())
    {
      CombinationSet_t combinations = graph.GetNodeCombinations(n, RecoType_t::kShower, 0, 22);
      for (auto const& combination : combinations)
      {
        if (_verbose)
        {
          std::cout << "[" << __FUNCTION__ << "] Comparing showers";
          for (auto const& particle : combination)
            std::cout << " " << graph.GetParticle(particle).ID();
        }
        double l = Likelihood(graph, combination);

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
    
    // output some info on particles found
    
    if (_verbose)
    {
      std::cout << "[" << __FUNCTION__ << "] " << event.size() << " objects found in event." << std::endl;;
      for (auto const& object : event)
      {
        std::cout << "[" << __FUNCTION__ << "] Gamma combination";
        for (auto const& id : object.first)
        {
          std::cout << " " << id;
        }
        std::cout << " with score " << object.second << std::endl;
      }
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
    if(_useEnergyPdf == false && _useMomentumPdf == false && _useMassPdf == false && _useRadLenPdf == false)
    {
      std::cout << "Error! Trying to perform a likelihood calculation, but no PDFs have been instantiated!" << std::endl;
      return 0;
    }
    
    double l = 1;
    
    // run through all the pdfs, and contribute to the score if the pdf is used

    if (_useEnergyPdf)
    {
      double energy = _extEnergy;
      for (auto const& particle : combination)
        energy += graph.GetParticle(particle).Energy();
      _energyVar->setVal(energy);
      double e = _energyPdf->getVal(*_energyVar);
      double emax = 1 / _energyPdf->getNorm(*_energyVar);
      double p = e / emax;
      l *= p;
    }
    
    if (_useMomentumPdf)
    {
      geoalgo::Vector momentum = *_extMomentum;
      for (auto const& particle : combination)
        momentum += graph.GetParticle(particle).Momentum();
      _momentumVar->setVal(momentum.Length());
      double e = _momentumPdf->getVal(*_momentumVar);
      double emax = 1 / _momentumPdf->getNorm(*_momentumPdf);
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
      
      // calculate invariant mass
      auto particle1 = graph.GetParticle(combination[0]);
      auto particle2 = graph.GetParticle(combination[1]);
      double angle = particle1.Momentum().Angle(particle2.Momentum());
      std::cout << "[" << __FUNCTION__ << "] DEBUG OUTPUT: ANGLE IS " << angle << std::endl;
    }
    
    return l;
  }
  
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
      
      if (_verbose)
        std::cout << "[" << __FUNCTION__ << "] Energy " << energy << " added to data set" << std::endl;
      
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
    if (_verbose)
      std::cout << "[" << __FUNCTION__ << "] Momentum " << momentum.Length() << "MeV added to data set" << std::endl;
    
    // tell the code you're using this pdf
    _useMomentumPdf = true;
    
    // set roofit variable & add to data set
    _momentumVar->setVal(momentum.Length());
    _momentumData->add(RooArgSet(*_momentumVar));
  }
  
  void NGammaBase::FillMassPdf(double mass)
  {
    if (_verbose)
      std::cout << "[" << __FUNCTION__ << "] Mass " << mass << "MeV added to data set" << std::endl;
    
    // tell the code you're using this pdf
    _useMassPdf = true;
    
    // set roofit variable & add to data set
    _massVar->setVal(mass);
    _massData->add(RooArgSet(*_massVar));
  }
  
  void NGammaBase::FillRadLenPdf(double rad_len)
  {
    // make sure the value is non-negative
    if (rad_len >= 0)
    {
      // tell the code you're using this pdf
      _useRadLenPdf = true;
      
      if (_verbose)
        std::cout << "[" << __FUNCTION__ << "] Particle radiation length " << rad_len << " added to data set" << std::endl;
      
      // set val & add to data set
      _radLenVar->setVal(rad_len);
      _radLenData->add(RooArgSet(*_radLenVar));
    }
    
    else
      if (_verbose)
        std::cout << "[" << __FUNCTION__ << "] Radiation length value " << rad_len << " cannot be negative!" << std::endl;
  }
  
  void NGammaBase::SaveParams()
  {
    // make sure at least one pdf has been instantiated before opening output file
    if (!(_useEnergyPdf || _useMomentumPdf || _useMassPdf || _useRadLenPdf))
    {
      std::cout << "[" << __FUNCTION__ << "] No PDFs were instantiated! Not saving output." << std::endl;
      return;
    }
    
    // save energy params
    if (_useEnergyPdf)
    {
      if (_verbose)
      {
        std::cout << "[" << __FUNCTION__ << "] ";
        _energyData->Print();
      }
      
      // get pointers to mean & sigma
      RooRealVar *meanVar, *sigmaVar;
      meanVar = (RooRealVar*)(_energyPdf->getVariables()->find("energyPdf_Gaus_mean"));
      sigmaVar = (RooRealVar*)(_energyPdf->getVariables()->find("energyPdf_Gaus_sigma"));
      
      // stop it from crashing somehow
      double min, max, mean;
      _energyData->getRange(*_energyVar, min, max);
      std::cout << "[" << __FUNCTION__ << "] MEAN VALUE IS " << _energyData->mean(*_energyVar) << std::endl;
      min = 0.9 * min;
      max = 1.1 * max;
      mean = _energyData->mean(*_energyVar);
      meanVar->setRange(0.1, max);
      sigmaVar->setRange(0.1, max);
      meanVar->setVal(mean);
      sigmaVar->setVal(mean);
      
      // fit pdf to data
      RooFitResult* energyResult = _energyPdf->fitTo(*_energyData, RooFit::Save(), RooFit::PrintLevel(-1));
      energyResult->Print();
      
      // plot pdf & data
      TCanvas *c = new TCanvas("c","",1000,500);
      RooPlot* frame = _energyVar->frame();
      _energyData->plotOn(frame, RooFit::MarkerColor(kRed), RooFit::LineColor(kRed));
      _energyPdf->plotOn(frame, RooFit::LineColor(kRed));
      frame->Draw();
      c->SaveAs("./energy_pdf.png");

      // save the parameter set to a config file (doesn't work right now)
      std::vector<double> darray;
      auto& params = OutputPSet();
      
      darray.push_back(meanVar->getVal());
      darray.push_back(sigmaVar->getVal());
      
      params.add_value("energy_params",::fcllite::VecToString(darray));
    }
    
    // save mass params
    if (_useMassPdf)
    {
      if (_verbose)
      {
        std::cout << "[" << __FUNCTION__ << "] ";
        _energyData->Print();
      }
      
      // get pointers to mean & sigma
      RooRealVar *meanVar, *sigmaVar;
      meanVar = (RooRealVar*)(_massPdf->getVariables()->find("massPdf_Gaus_mean"));
      sigmaVar = (RooRealVar*)(_massPdf->getVariables()->find("massPdf_Gaus_sigma"));
      
      // stop it from crashing somehow
      double min, max, mean;
      _massData->getRange(*_massVar, min, max);
      min *= 0.9;
      max *= 1.1;
      mean = _massData->mean(*_massVar);
      _massVar->setRange(0.1, 250);
      meanVar->setRange(0.1, max);
      sigmaVar->setRange(0.1, max);
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
      frame->GetXaxis()->SetLimits(0,200);
      frame->Draw();
      c->SaveAs("./mass_pdf.png");
      
      // save the parameter set to a config file (doesn't work right now)
      std::vector<double> darray;
      auto& params = OutputPSet();
      
      darray.push_back(meanVar->getVal());
      darray.push_back(sigmaVar->getVal());
      
      params.add_value("mass_params",::fcllite::VecToString(darray));
    }
  }
}

#endif
