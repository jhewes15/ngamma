import ROOT
from ROOT import ertool
from random import randint

outfile = ROOT.TFile("out.root","RECREATE")

mgr = ertool.Manager()
mgr.AddCfgFile('ngamma.cfg')

algo_empart = ertool.AlgoEMPart()
#algo_empart.setVerbose(True)

algo_pi0 = ertool.ERAlgoJeremyPi0()
algo_pi0.SetVerbose(True)
algo_pi0.TrainingMode(True)

mgr.AddAlgo(algo_empart)
mgr.AddAlgo(algo_pi0)

mgr.Initialize()

rand = ROOT.TRandom()

for i in xrange(100):
  mgr.ClearData()
  for j in xrange(randint(1,5)):
    shower = algo_pi0.Generate()
    myid = ertool.RecoInputID_t()
    mgr.Add(shower,myid)

  mgr.Process()

mgr.Finalize(outfile)
