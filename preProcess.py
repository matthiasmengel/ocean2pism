#!/usr/bin/env python

""" This script produces input for the runDiffuse.py script.
    it concatenates the source nc files, chooses the sigma levels and
    sets time axis and missing values. """

import sys, os, imp
import time
import PreAndPostProcess; reload(PreAndPostProcess)

sourcepath = "/iplex/01/tumble/mengel/pismSourceData/20120705_BriosHadgem2HistoricalRcp45Rcp85/origdata_rcp85/"
sourcepath = "/iplex/01/tumble/mengel/pismSourceData/20120705_BriosHadgem2HistoricalRcp45Rcp85/origdata_hist/"
#workpath   = "/iplex/01/tumble/mengel/pismInputData/"
workpath   = "/scratch/01/mengel/pismInputData/"
runname    = os.path.basename(os.getcwd())
outpath    = workpath + runname
## if files were already concatenated earlier
concat_done=True

#briosids = ["BRIOS.HadGem2_RCP85"]
briosids = ["BRIOS.HadGem2_20C"]
levelids = ["1"]#,"1,2","1,2,3","2","3","2,3"]
runids = []

for briosid in briosids:
  for levelid in levelids:
    #runids.append(briosid + "__" + levelid)

    #briosid, levelid = runid.split("__")
    #### choose sigma layer
    pp = PreAndPostProcess.PreAndPostProcess(sourcepath, outpath, briosid)
    if not concat_done:
      pp.concatenate()
    #pp.yearly_ave_taxis_missval_combine()
    pp.choose_sigmalevel(levelid) # e.g. "1,2,3" or "2"
