import ROOT
from ROOT import ertool
from random import randint

outfile = ROOT.TFile("out.root","RECREATE")

mgr = ertool.Manager()
mgr.AddCfgFile('ngamma.cfg')

algo_empart = ertool.AlgoEMPart()
#algo_empart.setVerbose(True)

algo_ngamma = ertool.ERAlgoNGamma()
#algo_ngamma.SetVerbose(True)

mgr.AddAlgo(algo_empart)
mgr.AddAlgo(algo_ngamma)

mgr.Initialize()

rand = ROOT.TRandom()

for i in xrange(100):
  mgr.ClearData()
  for j in xrange(randint(1,5)):
    shower = algo_ngamma.Generate()
    myid = ertool.RecoInputID_t()
    mgr.Add(shower,myid)

  mgr.Process()

mgr.Finalize(outfile)
