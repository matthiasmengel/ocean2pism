#!/usr/bin/env python
"""
This program allows temperature and salinity diffuse from points where it is defined
to regions where there is ice. Note that only regions below sea level are included
in calculation. matthias.mengel@pik

"""
#import matplotlib.pylab as plt
import os
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import time

scen  = "rcp85"
infile  = "../projected_files/brios_15km_" + scen + ".nc"
outfile = "../diffused_files/brios_15km_" + scen + ".nc"

pismfile = "/iplex/01/tumble/mengel/pismOut/pismDev_gh002_071ECEarthBoundsEsia56NoMass15km/NoMass.nc"

a         = 1e6          # Diffusion constant.
timesteps = 150000 # Number of time-steps to evolve system.

nci = nc.Dataset(infile,  'r')
ncp = nc.Dataset(pismfile,'r')
x = nci.variables['x'][:]
y = nci.variables['y'][:]
dx = x[1]-x[0]
dy = y[1]-y[0]

dx2=dx**2 # To save CPU cycles, we'll compute Delta x^2
dy2=dy**2 # and Delta y^2 only once and store them.

# For stability, this is the largest interval possible
# for the size of the time-step:
dt = dx2*dy2/( 2*a*(dx2+dy2) )

def evolve_ts2(ui, uii):
   """
   This function uses a numpy expression to
   evaluate the derivatives in the Laplacian, and
   calculates u[i,j] based on ui[i,j].
   """
   # diffusion
   u[1:-1, 1:-1] = ui[1:-1, 1:-1] + a*dt*( (ui[2:, 1:-1] - 
    2*ui[1:-1, 1:-1] + ui[:-2, 1:-1])/dx2 + (ui[1:-1, 2:] - 
    2*ui[1:-1, 1:-1] + ui[1:-1, :-2])/dy2 )
   # set known brios values back to initial
   u[notmask] = uii[notmask]
   # set values with topg>0 that are adjacent topg<0 cells
   # to the mean value of these topg<0 neighbours
   ud[1:-1, 1:-1] = ((ui[ :-2,  1:-1]*abovesea1[1:-1, 1:-1] + 
                      ui[2:,   1:-1] *abovesea2[1:-1, 1:-1] +
                      ui[1:-1,  :-2] *abovesea3[1:-1, 1:-1] +
                      ui[1:-1, 2:]   *abovesea4[1:-1, 1:-1]) / 
                     neighboursbelow[1:-1, 1:-1] )
   u[neighboursbelowmask] = ud[neighboursbelowmask]
   return u

belowsea = ncp.variables['topg'][0,:,:] <= 0.
abovesea = ncp.variables['topg'][0,:,:] >  0.
temp_all = nci.variables['thetao'][:,:,:]
salt_all = nci.variables['salinity'][:,:,:]
# these points are above sea and have a neighbour below
# we need these for a good boundary condition at margins
# to above sea level. above sea level temperatures should not 
# influence these below.
abovesea1 = np.copy(belowsea)
abovesea1[:] = False
abovesea1[1:-1,:] = (belowsea[:-2,:]) &  (abovesea[1:-1,:])
abovesea2 = np.copy(belowsea)
abovesea2[:] = False
abovesea2[1:-1,:] = (belowsea[2:,:]) &  (abovesea[1:-1,:])
abovesea3 = np.copy(belowsea)
abovesea3[:] = False
abovesea3[:,1:-1] = (belowsea[:,:-2]) &  (abovesea[:,1:-1])
abovesea4 = np.copy(belowsea)
abovesea4[:] = False
abovesea4[:,1:-1] = (belowsea[:,2:]) &  (abovesea[:,1:-1])

neighboursbelowmask = abovesea1 +abovesea2 +abovesea3 +abovesea4 
neighboursbelow     = abovesea1*1 +abovesea2*1 +abovesea3*1 +abovesea4*1 
innerpointsabovesl  = (neighboursbelowmask == 0) & (abovesea)

temp_out=np.copy(temp_all)
temp_out[:]=0.
temp_outraw=np.copy(temp_all)
temp_outraw[:]=0.
salt_out=np.copy(salt_all)
salt_out[:]=0.

print temp_out.nbytes/1000
print salt_out.nbytes/1000

f = open('./logfile', 'w')
f.write("start diffusion\n")
print "start diffusion\n"
for t in np.arange(0,temp_all.shape[0]):
    print "timestep" + str(t)+ "\n"
    f.write("timestep" + str(t)+ "\n") 
    temp_available = temp_all[t,:,:]
    salt_available = salt_all[t,:,:]
    # temperature
    uii = ma.copy(temp_available)
    notmask = ~uii.mask
    # get rid of mask
    ui  = np.copy(temp_available)
    ui[ui>500.]  = 273.15 -2.
    ui[ui<-500.] = 273.15 -2.
    u    = np.copy(ui)
    ud   = np.zeros(ui.shape)

    #thetao1   = np.zeros([timesteps,u.shape[0],u.shape[1]])
    #salinity1 = np.zeros([timesteps,u.shape[0],u.shape[1]])
    # Now, start the time evolution calculation...
    tstart = time.time()
    m=0
    while m < timesteps:
      u  = evolve_ts2(ui, uii)
      ui = u
    #  thetao1[m,:,:] = u
      #print "Computing u for m =" + str(m) + "\n"
      #f.write("Computing u for m =" + str(m) + "\n")
      m += 1
    tfinish = time.time()
    oceantemp = np.copy(u)
    #oceantemp[innerpointsabovesl] = 0.
    oceantemp_raw = np.copy(uii)

    # salinity
    uii = ma.copy(salt_available)
    notmask = ~uii.mask
    # get rid of mask
    ui  = np.copy(salt_available)
    ui[ui>500.] = 34.8
    ui[ui<-500.] = 34.8
    u    = np.copy(ui)
    uu   = np.copy(ui)
    ud   = np.zeros(ui.shape)

    # Now, start the time evolution calculation...
    tstart = time.time()
    f.write("start salinity\n")
    m=0
    while m < timesteps:
      u  = evolve_ts2(ui, uii)
      ui = u
      m += 1
    #  salinity1[m,:,:] = u
#      print "Computing u for m =", m
    tfinish = time.time()
    salinity = np.copy(u)
    temp_out[t,:,:] = oceantemp
    temp_outraw[t,:,:] = oceantemp_raw
    salt_out[t,:,:] = salinity
f.close()

ncout = nc.Dataset(outfile, 'w', format='NETCDF3_CLASSIC')
ncout.createDimension('time',size=None)
ncout.createDimension('x',size=len(nci.variables['x'][:]))
ncout.createDimension('y',size=len(nci.variables['y'][:]))
ncvart  = ncout.createVariable( 'thetao','float32',('time','y','x') )
ncvars = ncout.createVariable( 'salinity','float32',('time','y','x')  )
#ncvart  = ncout.createVariable( 'thetao','float32',('y','x') )
#ncvars = ncout.createVariable( 'salinity','float32',('y','x')  )
ncvartr = ncout.createVariable( 'thetao_raw','float32',('time','y','x')  )
nct   = ncout.createVariable( 'time','float32',('time',) )
ncx   = ncout.createVariable( 'x','float32',('x',) )
ncy   = ncout.createVariable( 'y','float32',('y',) )
nclon = ncout.createVariable( 'lon','float32',('y','x') )
nclat = ncout.createVariable( 'lat','float32',('y','x') )
ncthk = ncout.createVariable( 'thk','float32',('y','x') )
nctg  = ncout.createVariable( 'topg','float32',('y','x') )

#nct[:] = 0.
nct[:] = nci.variables['time'][:]
ncvart[:]  = temp_out
ncvars[:]  = salt_out
ncvartr[:] = temp_outraw
#ncvart[:]  = thetao1 + 273.15
#ncvars[:]  = salinity1

ncx[:]     = x
ncy[:]     = y
#ncvars[:]  = uu
ncthk[:]    = ncp.variables['thk'][:]
nctg[:]     = ncp.variables['topg'][:]

ncy.units = 'meters'
ncx.units = 'meters'
ncvart.units = 'Kelvin'
ncvartr.units = 'Kelvin'
ncvars.units = 'psu'
ncthk.units  = 'meters'
nctg.units   = 'meters'

ncout.close()
nci.close()
ncp.close()


