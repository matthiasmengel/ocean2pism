
import os
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import time

class DiffuseOcean:
  """A class to diffuse Ocean Data over land where topg < 0.
  """
  def __init__(self, infile, outfile, pismfile, timesteps):
    self.infile    = infile
    self.outfile   = outfile
    self.pismfile  = pismfile
    self.timesteps = timesteps
    self.a         = 1e6    # Diffusion constant.

#scen  = "rcp85"
#infile  = "../projected_files/brios_15km_" + scen + ".nc"
#outfile = "../diffused_files/brios_15km_" + scen + ".nc"
#pismfile = "/iplex/01/tumble/mengel/pismOut/pismDev_gh002_071ECEarthBoundsEsia56NoMass15km/NoMass.nc"
#timesteps = 150000 # Number of time-steps to evolve system.

  def getInputData(self):
    nci = nc.Dataset(self.infile,  'r')
    ncp = nc.Dataset(self.pismfile,'r')
    self.topg = ncp.variables['topg'][0,:,:]
    self.thk  = ncp.variables['thk'][0,:,:]
    self.temp_all = nci.variables['thetao'][:,:,:]
    self.salt_all = nci.variables['salinity'][:,:,:]
    self.time = nci.variables['time'][:]
    self.x = nci.variables['x'][:]
    self.y = nci.variables['y'][:]
    dx = self.x[1]-self.x[0]
    dy = self.y[1]-self.y[0]
    self.dx2=dx**2
    self.dy2=dy**2
    # For stability, this is the largest interval possible
    # for the size of the time-step:
    self.dt = self.dx2*self.dy2/( 2*self.a*(self.dx2+self.dy2) )
    nci.close()
    ncp.close()

  def findAboveAndBelowSea(self):
    belowsea = self.topg <= 0.
    abovesea = self.topg >  0.
    abovesea1 = np.copy(abovesea)
    abovesea1[:] = False
    abovesea2 = np.copy(abovesea1)
    abovesea3 = np.copy(abovesea1)
    abovesea4 = np.copy(abovesea1)
    abovesea1[1:-1,:] = (belowsea[:-2,:]) &  (abovesea[1:-1,:])
    abovesea2 = np.copy(abovesea1)
    abovesea2[1:-1,:] = (belowsea[2:,:]) &  (abovesea[1:-1,:])
    abovesea3[:,1:-1] = (belowsea[:,:-2]) &  (abovesea[:,1:-1])
    abovesea4[:,1:-1] = (belowsea[:,2:]) &  (abovesea[:,1:-1])
    self.neighboursbelowmask = abovesea1 +abovesea2 +abovesea3 +abovesea4
    self.neighboursbelow     = abovesea1*1 +abovesea2*1 +abovesea3*1 +abovesea4*1
    self.innerpointsabovesl  = (self.neighboursbelowmask == 0) & (abovesea)
    self.abovesea1 = abovesea1
    self.abovesea2 = abovesea2
    self.abovesea3 = abovesea3
    self.abovesea4 = abovesea4

  def runDiffusion(self):
  
    def evolve_ts2(ui, uii):
      """ This function uses a numpy expression to evaluate the derivatives
          in the Laplacian, and calculates u[i,j] based on ui[i,j]. """
      # diffusion
      u[1:-1, 1:-1] = ui[1:-1, 1:-1] + self.a*self.dt*( (ui[2:, 1:-1] -
        2*ui[1:-1, 1:-1] + ui[:-2, 1:-1])/self.dx2 + (ui[1:-1, 2:] -
        2*ui[1:-1, 1:-1] + ui[1:-1, :-2])/self.dy2 )
      # set known brios values back to initial
      u[notmask] = uii[notmask]
      # set values with topg>0 that are adjacent topg<0 cells
      # to the mean value of these topg<0 neighbours
      ud[1:-1, 1:-1] = ((ui[ :-2,  1:-1] *self.abovesea1[1:-1, 1:-1] +
                          ui[2:,   1:-1] *self.abovesea2[1:-1, 1:-1] +
                          ui[1:-1,  :-2] *self.abovesea3[1:-1, 1:-1] +
                          ui[1:-1, 2:]   *self.abovesea4[1:-1, 1:-1]) /
                        self.neighboursbelow[1:-1, 1:-1] )
      u[self.neighboursbelowmask] = ud[self.neighboursbelowmask]
      return u

    self.temp_out    = np.copy(self.temp_all)
    self.temp_out[:] = 0.   
    #temp_outraw = np.copy(temp_out)
    self.salt_out    = np.copy(self.temp_out)

    print "start diffusion\n"
    for t in np.arange(0,self.temp_all.shape[0]):
      print "timestep" + str(t)+ "\n"
      uii = ma.copy(self.temp_all[t,:,:])
      ui  = np.copy(self.temp_all[t,:,:]) # get rid of mask
      notmask = ~uii.mask
      ui[ui>500.]  = 273.15 -2.
      ui[ui<-500.] = 273.15 -2.
      u   = np.copy(ui)
      ud  = np.zeros(ui.shape)
      m=0
      while m < self.timesteps:
        print "Computing temp for m =" + str(m) + "\n"
        u  = evolve_ts2(ui, uii)
        ui = u
        m += 1
      self.temp_out[t,:,:] = u
      #temp_outraw[t,:,:] = uii

      uii = ma.copy(self.salt_all[t,:,:])
      ui  = np.copy(self.salt_all[t,:,:])  # get rid of mask
      notmask = ~uii.mask
      ui[ui>500.] = 34.8
      ui[ui<-500.] = 34.8
      u    = np.copy(ui)
      ud   = np.zeros(ui.shape)
      print "start salinity"
      m=0
      while m < self.timesteps:
        print "Computing salt for m =" + str(m) + "\n"
        u  = evolve_ts2(ui, uii)
        ui = u
        m += 1
      self.salt_out[t,:,:] = u

  def writeNetcdf(self):
    ncout = nc.Dataset(self.outfile, 'w', format='NETCDF3_CLASSIC')
    ncout.createDimension('time',size=None)
    ncout.createDimension('x',size=len(self.x))
    ncout.createDimension('y',size=len(self.y))
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

    nct[:]     = self.time
    ncvart[:]  = self.temp_out
    ncvars[:]  = self.salt_out
    ncvartr[:] = self.temp_all
    ncx[:]     = self.x
    ncy[:]     = self.y
    ncthk[:]    = self.thk
    nctg[:]     = self.topg

    ncy.units = 'meters'
    ncx.units = 'meters'
    ncvart.units = 'Kelvin'
    ncvartr.units = 'Kelvin'
    ncvars.units = 'psu'
    ncthk.units  = 'meters'
    nctg.units   = 'meters'

    ncout.close()
        