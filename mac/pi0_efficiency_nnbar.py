import ROOT
from ROOT import larlite as fmwk
from ROOT import ertool
from random import randint

for j in range(1,11):
  
  my_proc = fmwk.ana_processor()
  my_proc.set_io_mode(fmwk.storage_manager.kREAD)

  outfile = ROOT.TFile("out.root","RECREATE")

  my_proc.add_input_file('/Users/jhewes15/neutrino/larlite/data/nnbar/larlite_mcinfo.root')

  anaunit = fmwk.ExampleERSelection()

  anaunit._mgr.ClearCfgFile()
  anaunit._mgr.AddCfgFile("ngamma.cfg")

  anaunit._mgr._mc_for_ana = True

  algo_empart = ertool.AlgoEMPart()
  algo_empart.setVerbose(False)

  algo_pi0 = ertool.ERAlgoJeremyPi0()
  algo_pi0.SetVerbose(False)
  algo_pi0.TrainingMode(False)

  ana_pi0 = ertool.ERAnaJeremyPi0()

  anaunit._mgr.AddAlgo(algo_empart)
  anaunit._mgr.AddAlgo(algo_pi0)
  anaunit._mgr.AddAna(ana_pi0)

  anaunit.SetShowerProducer(True,"mcreco")

  my_proc.add_process(anaunit)

  algo_pi0.SetCut(j)
  ana_pi0.SetBin(j)
  my_proc.run(j)

