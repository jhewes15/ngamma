import ROOT
from ROOT import larlite as fmwk
from ROOT import ertool
from random import randint

my_proc = fmwk.ana_processor()
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

my_proc.add_input_file('/Users/jhewes15/neutrino/larlite/data/nnbar/larlite_mcinfo.root')

anaunit = fmwk.ExampleERSelection()

algo_empart = ertool.AlgoEMPart()
#algo_empart.setVerbose(True)

algo_nnbar = ertool.ERAlgoNNbar()
algo_nnbar.SetVerbose(True)
algo_nnbar.TrainingMode(True)

ana_nnbar = ertool.ERAnaNNbar()

anaunit._mgr.AddAlgo(algo_empart)
anaunit._mgr.AddAlgo(algo_nnbar)
#anaunit._mgr.AddAna(ana_nnbar)

anaunit.SetShowerProducer(True,"mcreco")

my_proc.add_process(anaunit)
my_proc.run()

anaunit._mgr.StorePSet("nnbar.cfg")

