import sys, glob, os
import subprocess
sys.path.append("/home/mengel/gitProjects/DiffuseAntarcticOcean")

#import PreAndPostprocess; reload(PreAndPostprocess)

def concatenate_brios(inpath, infilepattern, outpath, addstring=""):
  outp    = outpath + "/concatenated/"
  infiles = glob.glob(inpath + "/" + infilepattern)
  outfile = infilepattern.replace(".*","")
  os.system("mkdir -p " + outp)

  infilestr = ""
  for fl in infiles:
    infilestr += fl + " "

  cmd = "ncrcat -O " + infilestr + " " + outp + outfile
  print "## create " + outp + outfile
  subprocess.check_call(cmd, shell=True)


def yearly_ave_taxis_missval(infile, outpath, addstring=""):
  outp    = outpath + "/yearmeantaxisandmissing/"
  starty  = "1950" if "20C" in infile else "2000"
  missval = "35" if ".s." in infile else "0"

  os.system("mkdir -p " + outp)

  cmd = "cdo yearmean -settaxis," + starty + "-01-01,00:00:00,1mon -setctomiss," + missval
  cmd += " " + outpath + "/concatenated/" + infile + " " + outp + infile

  print "## create " + outp + infile
  subprocess.check_call(cmd, shell=True)


def combine_fields(tfile, sfile, outpath, addstring=""):
  outp    = outpath + "/combinedfields/"

  os.system("mkdir -p " + outp)

  tf = outpath + "/yearmeantaxisandmissing/" + tfile
  sf = outpath + "/yearmeantaxisandmissing/" + sfile
  outf = tfile.replace(".t.",".")
  cmd = "cdo -O merge " + tf + " " + sf + " " + outp + outf

  print "## create " + outp + outf
  subprocess.check_call(cmd, shell=True)


def choose_sigmalevel_rename(infile, outpath, level, inpath=""):
  outp    = outpath + "/sigmalevel/"

  inpath = outpath if inpath=="" else inpath

  os.system("mkdir -p " + outp)
  infl = inpath + "/combinedfields/" + infile
  cmd = "cdo -O -vertmean -sellevel," + str(level) + " -chname,temperature,thetao " + infl + " " + outp + infile

  print "## create " + outp + infile
  subprocess.check_call(cmd, shell=True)

  cmd = "ncatted -O -a units,salinity,o,c,g/kg " + outp + infile
  subprocess.check_call(cmd, shell=True)

def average_over_years(infile, outpath, tstep1, tsteplast, addstring=""):
  outp    = outpath + "/multiyearaverages/"

  os.system("mkdir -p " + outp)
  infl = outpath + "/diffused/" + infile
  cmd = "ncra -F -O -d time," + str(tstep1) + "," + str(tsteplast) + ",1 " + infl + " " + outp + infile
  print "## create " + outp + infile
  subprocess.check_call(cmd, shell=True)
  # this has to be in the diffusion code soon
  cmd = "ncatted -O -a units,salinity,o,c,g/kg " + outp + infile
  subprocess.check_call(cmd, shell=True)

  #starty  = "1950" if "20C" in infile else "2000"

  #cmd = "cdo yearmean -settaxis," + starty + "-01-01,00:00:00,1mon
  #cmd += " " + outpath + "/concatenated/" + infile + " " + outp + infile

  #print "## create " + outp + infile
  #subprocess.check_call(cmd, shell=True)



