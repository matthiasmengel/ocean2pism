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
diffuse_variables = ["thetao", "salinity", "ismelt"]
diffuse_missvals  = {"thetao":-2., "salinity":34.8, "ismelt":0.}
#workpath   = "/iplex/01/tumble/mengel/pismInputData/"
workpath   = "/scratch/01/mengel/pismInputData/"
destination_grid_file = "/iplex/01/tumble/mengel/pismOut/pismDev_gh002_071ECEarthBoundsEsia56NoMass15km/NoMass.nc"
runname    = os.path.basename(os.getcwd())
outpath    = workpath + runname
infile     = outpath + "/sigmalevel/" + "BRIOS.HadGem2_20C_sigma1.nc"
outfile    = outpath + "/" + os.path.basename(infile).strip(".nc") + "_diffused" + str(diffuse_timesteps) + ".nc"

# only copy data necessary for pism if lite
lite = True

os.system("mkdir -p " + outpath + "/diffused/")
print "initialize"
dd = DiffuseOcean.DiffuseOcean(infile, outfile, destination_grid_file, diffuse_timesteps, diffuse_variables, diffuse_missvals)
dd.getTimeIndependentInputData()
dd.prepareProjection()
dd.findAboveAndBelowSea()
dd.findAboveAndBelowSea()


# try to get the communicator object to see whether mpi is available:
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    available = False if comm.size < 2 else True
    print "mpi available, with ", comm.size, " processors."
except:
    available = False

if available:
  result_parallel = mpi4py_map.map(dd.projAndDiffu, xrange(dd.nctimesteps), debug=True)
  dd.writeNetcdf(result_parallel, lite)
else:
  result_serial = map(dd.projAndDiffu, xrange(dd.nctimesteps))
  dd.writeNetcdf(result_serial, lite)

