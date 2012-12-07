
import os, time
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import interp
from pyproj import Proj

class DiffuseOcean:
  """A class to diffuse Ocean Data over land where topg < 0.
  """
  def __init__(self, infile, outfile, pismfile, timesteps):
    self.infile    = infile
    self.outfile   = outfile
    self.pismfile  = pismfile
    self.timesteps = timesteps
    self.a         = 1e6    # Diffusion constant.

  def getInputData(self):
    nci = nc.Dataset(self.infile,  'r')
    ncp = nc.Dataset(self.pismfile,'r')
    self.topg = ncp.variables['topg'][0,:,:].T
    self.thk  = ncp.variables['thk'][0,:,:].T
    self.olat = nci.variables['lat'][:,0]
    # take out last 3 longitude values, its double data
    self.olon = nci.variables['lon'][0,0:-3]
    self.otemp = nci.variables["thetao"][:,:,0:-3]
    self.osalt = nci.variables["salinity"][:,:,0:-3]
    self.time = nci.variables['time'][:]
    self.x = ncp.variables['x'][:]
    self.y = ncp.variables['y'][:]
    dx = self.x[1]-self.x[0]
    dy = self.y[1]-self.y[0]
    self.dx2=dx**2
    self.dy2=dy**2
    # For stability, this is the largest interval possible
    # for the size of the time-step:
    self.dt = self.dx2*self.dy2/( 2*self.a*(self.dx2+self.dy2) )
    nci.close()
    ncp.close()

  def projectOnPismGrid(self):

    def extend_interp(datafield):
      dfield_ext = ma.concatenate([ma.column_stack(southernlimitmask), datafield], 0)
      return interp(dfield_ext, self.olon, olat_ext, self.pismlon, self.pismlat)

    olat_ext = np.append(-82.1,self.olat)
    # add masked values at southernmost end
    southernlimitmask = ma.masked_all(len(self.olon))
    xgrid, ygrid = np.meshgrid(self.x,self.y)
    newprojection  = Proj(proj='stere',lat_0=-90,lon_0=0,lat_ts=-71,ellps='WGS84')
    # these should be the same as the lon lat variables in Le Brocq
    self.pismlon, self.pismlat = newprojection(xgrid,ygrid,inverse=True)
    self.projtemp = ma.zeros([len(self.time),xgrid.shape[0],xgrid.shape[1]])
    self.projsalt = ma.zeros([len(self.time),xgrid.shape[0],xgrid.shape[1]])

    for t in np.arange(0,len(self.time)):
      print "project timestep" + str(t)
      self.projtemp[t,:,:] = extend_interp(self.otemp[t,:,:])
      self.projsalt[t,:,:] = extend_interp(self.osalt[t,:,:])

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

    def run_diffuse(dfield, setmissval):

      def diffuse(ui, uii):
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

      uii = ma.copy(dfield)
      notmask = ~uii.mask
      self.notmask = notmask
      ui  = np.copy(uii) # get rid of mask
      ui[uii.mask] = setmissval
      u   = np.copy(ui)
      ud  = np.zeros(ui.shape)
      m=0
      while m < self.timesteps:
        print "Diffuse for m =" + str(m) + "\n"
        u  = diffuse(ui, uii)
        ui = u
        m += 1
      return u

    self.dfutemp    = ma.copy(self.projtemp)
    self.dfutemp[:] = 0.
    self.dfusalt    = ma.copy(self.dfutemp)

    print "start diffusion\n"
    for t in np.arange(0,self.projtemp.shape[0]):
      print "diffuse timestep " + str(t)+ "\n"
      self.dfutemp[t,:,:] = run_diffuse(self.projtemp[t,:,:],273.15 -2.0)
      self.dfusalt[t,:,:] = run_diffuse(self.projsalt[t,:,:],34.8)

  def writeNetcdf(self):
    ncout = nc.Dataset(self.outfile, 'w', format='NETCDF3_CLASSIC')
    ncout.createDimension('time',size=None)
    ncout.createDimension('x',size=len(self.x))
    ncout.createDimension('y',size=len(self.y))
    ncvart  = ncout.createVariable( 'thetao','float32',('time','y','x') )
    ncvars = ncout.createVariable( 'salinity','float32',('time','y','x')  )
    ncvartr = ncout.createVariable( 'thetao_raw','float32',('time','y','x')  )
    nct   = ncout.createVariable( 'time','float32',('time',) )
    ncx   = ncout.createVariable( 'x','float32',('x',) )
    ncy   = ncout.createVariable( 'y','float32',('y',) )
    nclon = ncout.createVariable( 'lon','float32',('y','x') )
    nclat = ncout.createVariable( 'lat','float32',('y','x') )
    ncthk = ncout.createVariable( 'thk','float32',('y','x') )
    nctg  = ncout.createVariable( 'topg','float32',('y','x') )
    ncx1  = ncout.createVariable( 'abovesea1','float32',('y','x') )
    ncx2  = ncout.createVariable( 'abovesea2','float32',('y','x') )
    ncx3  = ncout.createVariable( 'abovesea3','float32',('y','x') )
    ncx4  = ncout.createVariable( 'abovesea4','float32',('y','x') )
    ncx5  = ncout.createVariable( 'neighboursbelowmask','float32',('y','x') )
    ncx6  = ncout.createVariable( 'neighboursbelow','float32',('y','x') )
    ncx7  = ncout.createVariable( 'notmask','float32',('y','x') )
    nct[:]     = self.time
    ncvart[:]  = self.dfutemp
    ncvars[:]  = self.dfusalt
    ncvartr[:] = self.projtemp
    ncx[:]     = self.x
    ncy[:]     = self.y
    nclat[:]   = self.pismlat
    nclon[:]   = self.pismlon
    ncthk[:]    = self.thk
    nctg[:]     = self.topg
    ncx1[:]     = self.abovesea1
    ncx2[:]     = self.abovesea2
    ncx3[:]     = self.abovesea3
    ncx4[:]     = self.abovesea4
    ncx5[:]     = self.neighboursbelowmask
    ncx6[:]     = self.neighboursbelow
    ncx7[:]     = self.notmask
    ncy.units = 'meters'
    ncx.units = 'meters'
    ncvart.units = 'Kelvin'
    ncvartr.units = 'Kelvin'
    ncvars.units = 'psu'
    ncthk.units  = 'meters'
    nctg.units   = 'meters'
    nct.units    = "years since 01-01-00"
    nct.calendar = "365_day"
    ncout.comment = "diffused over " + str(self.timesteps) + "."
    ncout.close()

