#!/usr/bin/env python
"""
This program allows temperature and salinity diffuse from points where it is defined
to regions where there is ice. Note that only regions below sea level are included
in calculation. matthias.mengel@pik

"""
import sys, os, imp
import time
import mpi
PreAndPostProcess = imp.load_source('PreAndPostProcess', '/home/mengel/gitProjects/DiffuseAntarcticOcean/PreAndPostprocess.py')
DiffuseOcean      = imp.load_source('DiffuseOcean', '/home/mengel/gitProjects/DiffuseAntarcticOcean/DiffuseOcean.py')

sourcepath = "/iplex/01/tumble/mengel/pismSourceData/20120730_BriosHadGem2DrivenHistorical"
workpath   = "/iplex/01/tumble/mengel/pismInputData/"
pismfile   = "/iplex/01/tumble/mengel/pismOut/pismDev_gh002_071ECEarthBoundsEsia56NoMass15km/NoMass.nc"
runname    = os.path.basename(os.getcwd())
outpath    = workpath + runname

runids = ["BRIOS.HadGem2_20C_II","BRIOS.HadGem2_20C_III", "BRIOS.HadGem2_RCP26", "BRIOS.HadGem2_RCP45"]
diffuse_timesteps = 5

runid = runids[0]

def processdata(runid):
  #### prepare
  pp = PreAndPostProcess.PreAndPostProcess(sourcepath, outpath, runid)
  #pp.concatenate()
  #pp.yearly_ave_taxis_missval_combine()
  pp.choose_sigmalevel("2") # e.g. "1,2,3" or "2"

  #### diffuse
  infile  = outpath + "/sigmalevel/" + runid + ".nc"
  outfile = outpath + "/diffused/" + runid + ".nc"
  os.system("mkdir -p " + outpath + "/diffused/")
  print "initialize"
  dd = DiffuseOcean.DiffuseOcean(infile, outfile, pismfile, diffuse_timesteps)
  print "get input data"
  dd.getInputData()
  print "project on pism grid"
  dd.projectOnPismGrid()
  print "find points above sl"
  dd.findAboveAndBelowSea()
  print "run diffusion"
  dd.runDiffusion()
  print "write netcdf"
  dd.writeNetcdf()
  return 0


def master():
  for runid in runids:
    print runid
    result = {}
    mpi.submit_call("processdata", (runid,), id=runid)
  for runid in runids:
    result[runid] = mpi.get_result(runid)
  for runid in runids:
    print "------ "  + runid + ": ------\n" + str(result[runid])

x1 =time.strftime('%s')
mpi.run()
x2 =time.strftime('%s')
print 'it took ', int(x2)-int(x1),'seconds'
