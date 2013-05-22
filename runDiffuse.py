#!/usr/bin/env python
"""
This program allows temperature and salinity diffuse from points where it is defined
to regions where there is ice. Note that only regions below sea level are included
in calculation. matthias.mengel@pik

"""
import sys, os, imp
import time
import DiffuseOcean; reload(DiffuseOcean)
import numpy as np
import mpi4py_map

diffuse_timesteps = 1
#workpath   = "/iplex/01/tumble/mengel/pismInputData/"
workpath   = "/scratch/01/mengel/pismInputData/"
destination_grid_file = "/iplex/01/tumble/mengel/pismOut/pismDev_gh002_071ECEarthBoundsEsia56NoMass15km/NoMass.nc"
runname    = os.path.basename(os.getcwd())
outpath    = workpath + runname
infile     = outpath + "/sigmalevel/" + "BRIOS.HadGem2_20C__1.nc"
outfile    = outpath + "/" + os.path.basename(infile).strip(".nc") + "_diffused_" + str(diffuse_timesteps) + ".nc"

# only copy data necessary for pism if lite
lite = True

os.system("mkdir -p " + outpath + "/diffused/")
print "initialize"
dd = DiffuseOcean.DiffuseOcean(infile, outfile, destination_grid_file, diffuse_timesteps)
print "get input data"
#try:
dd.getInputData()
#except RuntimeError as error:
  #print "###" + runid + " failed, " + str(error)
  #sys.exit(0)
print "project on pism grid"
dd.projectOnPismGrid()
print "find points above sl"
dd.findAboveAndBelowSea()
print "run diffusion"
dd.runDiffusion()
print "write netcdf"
dd.writeNetcdf(lite)



#def processdata(runid, pp):
  #briosid, levelid = runid.split("__")
  ##### choose sigma layer
  #pp = PreAndPostProcess.PreAndPostProcess(sourcepath, outpath, briosid)
  #pp.choose_sigmalevel(levelid) # e.g. "1,2,3" or "2"

  ##### diffuse
  #filestr = briosid + "__" + levelid.replace(",","") + ".nc"
  #infile  = outpath + "/sigmalevel/" + filestr
  #outfile = outpath + "/diffused/" + filestr

  #return 0


#def master():

  #for briosid in briosids:
    #pp = PreAndPostProcess.PreAndPostProcess(sourcepath, outpath, briosid)
    #pp.concatenate()
    #pp.yearly_ave_taxis_missval_combine()

  #for runid in runids:
    #print runid
    #result = {}
    #mpi.submit_call("processdata", (runid,pp,), id=runid)
  #for runid in runids:
    #result[runid] = mpi.get_result(runid)
  #for runid in runids:
    #print "------ "  + runid + ": ------\n" + str(result[runid])

#x1 =time.strftime('%s')
#mpi.run()
#x2 =time.strftime('%s')
#print 'it took ', int(x2)-int(x1),'seconds'
