import sys, glob, os
import subprocess
sys.path.append("/home/mengel/gitProjects/DiffuseAntarcticOcean")

import PreAndPostprocess; reload(PreAndPostprocess)

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


def choose_sigmalevel(infile, outpath, level, addstring=""):
  outp    = outpath + "/sigmalevel/"

  os.system("mkdir -p " + outp)
  infl = outpath + "/combinedfields/" + infile
  cmd = "cdo -O sellevel," + str(level) + " " + infl + " " + outp + infile

  print "## create " + outp + infile
  subprocess.check_call(cmd, shell=True)


