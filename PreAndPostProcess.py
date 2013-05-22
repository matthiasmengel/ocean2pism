import sys, glob, os
import subprocess

class PreAndPostProcessException(Exception):
  pass

class PreAndPostProcess:

  def __init__(self, sourcepath, outpath, briosid):

    self.sourcepath = sourcepath
    self.outpath    = outpath
    self.briosid      = briosid
    self.concatp    = outpath + "/concatenated/"
    self.yearavep   = outpath + "/yearmeantaxisandmissing/"
    self.combinep   = outpath + "/combinedfields/"
    self.sigmalevp  = outpath + "/sigmalevel/"
    self.timeavep   = outpath + "/timeave/"

    self._makedirs()

  def _makedirs(self):

    for directory in [self.concatp, self.yearavep, self.combinep, self.sigmalevp]:
      os.system("mkdir -p " + directory + " 2> /dev/null")
      #os.system("rm " + directory + "/*.nc* 2> /dev/null")


  def concatenate(self):

    for pattern in [".t", ".s", ".ismelt"]:
      #print self.sourcepath + "/" + self.briosid + pattern + ".*.nc"
      infiles  = sorted(glob.glob(self.sourcepath + "/" + self.briosid + pattern + ".*.nc"))
      #print infiles
      if len(infiles) == 0:
        raise PreAndPostProcessException("no files to concatenate in " + self.sourcepath)

      concfile = self.concatp + "/" + self.briosid + pattern + ".nc"
      #print infiles
      infilestr = ""
      for fl in infiles:
        infilestr += fl + " "

      cmd = "ncrcat -O " + infilestr + " " + concfile
      print "## create " + concfile
      subprocess.check_call(cmd, shell=True)


  def yearly_ave_taxis_missval_combine(self):

    timeavefiles = ""
    for pt in [".t", ".s", ".ismelt"]:
      starty  = "1950" if "20C" in self.briosid else "2000"
      missval = "35" if pt == ".s" else "0"

      timeavefile = self.yearavep + self.briosid + pt + ".nc"
      cmd = "cdo yearmean -settaxis," + starty + "-01-01,00:00:00,1mon -setctomiss," + missval
      cmd += " " + self.concatp + self.briosid + pt + ".nc " + timeavefile

      print "## " + cmd
      subprocess.check_call(cmd, shell=True)
      timeavefiles += timeavefile + " "

    enlargefile = self.yearavep + self.briosid + ".ismelt.nc"

    mergefile = self.combinep + self.briosid + ".nc"
    cmd = "cdo -O merge " + timeavefiles + " " + mergefile
    print "## " + cmd
    subprocess.check_call(cmd, shell=True)


  def choose_sigmalevel(self, level, inputpath=""):

    inp = self.combinep if inputpath =="" else inputpath

    levstr = level.replace(",","")
    infl = inp  + self.briosid + ".nc"
    meanf = self.sigmalevp + self.briosid + "__" + levstr + "ts.nc"
    cmd = "cdo -O -vertmean -sellevel," + str(level) + " -chname,temperature,thetao " + infl + " " + meanf
    #cmd = "cdo -O -vertmean -chname,temperature,thetao " + infl + " " + meanf
    print "## " + cmd
    subprocess.check_call(cmd, shell=True)

    # change to pism understandable salinity units
    cmd = "ncatted -O -a units,salinity,o,c,g/kg " + meanf
    subprocess.check_call(cmd, shell=True)

    cmd = "ncatted -O -a extracted_sigma_level,global,a,c," + str(level) + " " + meanf
    subprocess.check_call(cmd, shell=True)

    # add again the ismelt as it was lost by sellevel
    outf = self.sigmalevp + self.briosid + "__" + levstr + ".nc"
    cmd = "cdo -O merge " + meanf + " " + self.yearavep + self.briosid + ".ismelt.nc " + outf
    subprocess.check_call(cmd, shell=True)


  def average_over_years(self, levelid, year1, yearend):
    outp   = self.timeavep
    levstr = levelid.replace(",","")
    infl   = self.outpath + "/diffused/" + self.briosid + "__" + levstr + ".nc"
    outf   = outp + self.briosid + "__" + levstr + "_" + str(year1) + "_" + str(yearend) + ".nc"
    os.system("mkdir -p " + outp)

    #cmd = "ncra -F -O -d time," + str(tstep1) + "," + str(tsteplast) + ",1 " + infl + " " + outf
    cmd = "cdo timmean -selyear," + str(year1) + "/" + str(yearend) + " " + infl + " " + outf
    print "## " + cmd
    subprocess.check_call(cmd, shell=True)
    ## this has to be in the diffusion code soon
    #cmd = "ncatted -O -a units,salinity,o,c,g/kg " + outp + infile
    #subprocess.check_call(cmd, shell=True)


  #def combine_fields(self, tfile, sfile, outpath, addstring=""):
    #outp    = outpath + "/combinedfields/"

    #os.system("mkdir -p " + outp)

    #tf = outpath + "/yearmeantaxisandmissing/" + tfile
    #sf = outpath + "/yearmeantaxisandmissing/" + sfile
    #outf = tfile.replace(".t.",".")
    #cmd = "cdo -O merge " + tf + " " + sf + " " + outp + outf

    #print "## create " + outp + outf
    #subprocess.check_call(cmd, shell=True)




