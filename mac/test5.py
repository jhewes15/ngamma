import ROOT
from ROOT import larlite as fmwk
from ROOT import ertool
from random import randint

my_proc = fmwk.ana_processor()
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

outfile = ROOT.TFile("out.root","RECREATE")

for i in range(1,201):
  my_proc.add_input_file('/Users/jhewes15/neutrino/larlite/data/pi0/larlite_mcinfo_{}.root'.format(i))

anaunit = fmwk.ExampleERSelection()

anaunit._mgr.ClearCfgFile()
anaunit._mgr.AddCfgFile("ngamma.cfg")

anaunit._mgr._mc_for_ana = True

algo_empart = ertool.AlgoEMPart()
#algo_empart.setVerbose(True)

algo_pi0 = ertool.ERAlgoJeremyPi0()
#algo_pi0.SetVerbose(True)
algo_pi0.TrainingMode(False)

ana_pi0 = ertool.ERAnaJeremyPi0()

anaunit._mgr.AddAlgo(algo_empart)
anaunit._mgr.AddAlgo(algo_pi0)
anaunit._mgr.AddAna(ana_pi0)

anaunit.SetShowerProducer(True,"mcreco")

my_proc.add_process(anaunit)
my_proc.run()

