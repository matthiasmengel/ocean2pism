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

sourcepath = "/iplex/01/tumble/mengel/pismSourceData/20120705_BriosHistoricalRcp45Rcp85/origdata_rcp85/"
#workpath   = "/iplex/01/tumble/mengel/pismInputData/"
workpath   = "/scratch/01/mengel/pismInputData/"
pismfile   = "/iplex/01/tumble/mengel/pismOut/pismDev_gh002_071ECEarthBoundsEsia56NoMass15km/NoMass.nc"
runname    = os.path.basename(os.getcwd())
outpath    = workpath + runname

briosids = ["BRIOS.HadGem2_RCP85"]
levelids = ["1"]#,"1,2","1,2,3","2","3","2,3"]
runids = []
for briosid in briosids:
  for levelid in levelids:
    runids.append(briosid + "__" + levelid)

#runids = ["BRIOS.HadGem2_RCP85"]
diffuse_timesteps = 1000

concat_done=False

def processdata(runid, pp):
  briosid, levelid = runid.split("__")
  #### choose sigma layer
  pp = PreAndPostProcess.PreAndPostProcess(sourcepath, outpath, briosid)
  pp.choose_sigmalevel(levelid) # e.g. "1,2,3" or "2"

  #### diffuse
  filestr = briosid + "__" + levelid.replace(",","") + ".nc"
  infile  = outpath + "/sigmalevel/" + filestr
  outfile = outpath + "/diffused/" + filestr
  os.system("mkdir -p " + outpath + "/diffused/")
  print "initialize"
  dd = DiffuseOcean.DiffuseOcean(infile, outfile, pismfile, diffuse_timesteps)
  print "get input data"
  try:
    dd.getInputData()
  except RuntimeError as error:
    print "###" + runid + " failed, " + str(error)
    return 1
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

  for briosid in briosids:
    pp = PreAndPostProcess.PreAndPostProcess(sourcepath, outpath, briosid)
    pp.concatenate()
    pp.yearly_ave_taxis_missval_combine()

  for runid in runids:
    print runid
    result = {}
    mpi.submit_call("processdata", (runid,pp,), id=runid)
  for runid in runids:
    result[runid] = mpi.get_result(runid)
  for runid in runids:
    print "------ "  + runid + ": ------\n" + str(result[runid])

x1 =time.strftime('%s')
mpi.run()
x2 =time.strftime('%s')
print 'it took ', int(x2)-int(x1),'seconds'
