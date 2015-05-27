#!/usr/bin/env python

# This file is part of ocean2pism.

# ocean2pism is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ocean2pism is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ocean2pism.  If not, see <http://www.gnu.org/licenses/>.

""" This script produces input for the runDiffuse.py script.
    it concatenates the source nc files, chooses the sigma levels and
    sets time axis and missing values. """

import sys, os, imp
import time
import PreAndPostProcess; reload(PreAndPostProcess)

# sourcepath = "/iplex/01/tumble/mengel/pismSourceData/20120705_BriosHadgem2HistoricalRcp45Rcp85/origdata_rcp85/"
sourcepath = "/iplex/01/tumble/mengel/pismSourceData/20120705_BriosHadgem2HistoricalRcp45Rcp85/origdata_hist/"
#workpath   = "/iplex/01/tumble/mengel/pismInputData/"
workpath   = "/scratch/01/mengel/pismInputData/"
runname    = os.path.basename(os.getcwd())
outpath    = workpath + runname
## if files were already concatenated earlier
concat_done=False

#briosids = ["BRIOS.HadGem2_RCP85"]
briosids = ["BRIOS.HadGem2_20C"]
levelids = ["24"]#,"1,2","1,2,3","2","3","2,3"]
runids = []

for briosid in briosids:
  for levelid in levelids:

    pp = PreAndPostProcess.PreAndPostProcess(sourcepath, outpath, briosid)
    if not concat_done:
      pp.concatenate()
    pp.combine_fields_year_ave()
    pp.choose_sigmalevel(levelid) # e.g. "1,2,3" or "2"
