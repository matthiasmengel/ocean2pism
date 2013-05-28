
import os, time
import netCDF4 as nc
import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import interp
from pyproj import Proj
import datetime

class DiffuseOcean:
  """A class to diffuse Ocean Data over land where topg < 0.
  """
  def __init__(self, infile, outfile, destination_grid_file, diffu_t, diffuse_vars, diffuse_missvals):

    print infile
    self.infile    = infile
    self.outfile   = outfile
    self.destination_grid_file  = destination_grid_file
    self.diffu_t = diffu_t
    self.diffuse_vars = diffuse_vars
    self.diffuse_missvals = diffuse_missvals
    self.a         = 1e6    # Diffusion constant.


  def getTimeIndependentInputData(self):

    nci = nc.Dataset(self.infile,  'r')
    ncp = nc.Dataset(self.destination_grid_file,'r')
    self.topg = ncp.variables['topg'][0,:,:]
    self.thk  = ncp.variables['thk'][0,:,:]
    self.olat = np.squeeze(nci.variables['lat'][:])[:,0]
    # take out last 3 longitude values, its double data in Brios
    self.olon = np.squeeze(nci.variables['lon'][:])[0,0:-3]
    #self.otemp = np.squeeze(nci.variables["thetao"][:])[:,:,0:-3] +273.15
    #self.osalt = np.squeeze(nci.variables["salinity"][:])[:,:,0:-3]
    #self.omelt = np.squeeze(nci.variables["ismelt"][:])[:,:,0:-3]
    self.time = nci.variables['time'][:]
    self.timeunits = nci.variables['time'].units
    self.calendar  = nci.variables['time'].calendar
    self.x = ncp.variables['x'][:]
    self.y = ncp.variables['y'][:]
    dx = self.x[1]-self.x[0]
    dy = self.y[1]-self.y[0]
    self.dx2=dx**2
    self.dy2=dy**2
    # For stability, this is the largest interval possible
    # for the size of the time-step:
    self.dt = self.dx2*self.dy2/( 2*self.a*(self.dx2+self.dy2) )
    self.nctimesteps = len(self.time)
    self.history = nci.history
    nci.close()
    ncp.close()


  def prepareProjection(self):

    xgrid, ygrid = np.meshgrid(self.x,self.y)
    newprojection  = Proj(proj='stere',lat_0=-90,lon_0=0,lat_ts=-71,ellps='WGS84')
    # these should be the same as the lon lat variables in Le Brocq
    self.pismlon, self.pismlat = newprojection(xgrid,ygrid,inverse=True)

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


  def projAndDiffu(self, tstep):

    print "tstep",tstep

    def extend_interp(datafield):
      # add masked values at southernmost end
      southernlimitmask = ma.masked_all(len(self.olon))
      olat_ext          = np.append(-82.1,self.olat)
      print southernlimitmask.shape, datafield.shape
      dfield_ext = ma.concatenate([ma.column_stack(southernlimitmask), datafield], 0)
      return interp(dfield_ext, self.olon, olat_ext, self.pismlon, self.pismlat)

    def run_diffuse(diffuse_var):
      #for diffuse_var in self.diffuse_vars:

      def diffuse(ui, projdata):
        """ This function uses a numpy expression to evaluate the derivatives
            in the Laplacian, and calculates u[i,j] based on ui[i,j]. """
        # diffusion
        u[1:-1, 1:-1] = ui[1:-1, 1:-1] + self.a*self.dt*( (ui[2:, 1:-1] -
          2*ui[1:-1, 1:-1] + ui[:-2, 1:-1])/self.dx2 + (ui[1:-1, 2:] -
          2*ui[1:-1, 1:-1] + ui[1:-1, :-2])/self.dy2 )
        # set known brios values back to initial
        u[notmask] = projdata[notmask]
        # set values with topg>0 that are adjacent topg<0 cells
        # to the mean value of these topg<0 neighbours
        ud[1:-1, 1:-1] = ((ui[ :-2,  1:-1] *self.abovesea1[1:-1, 1:-1] +
                            ui[2:,   1:-1] *self.abovesea2[1:-1, 1:-1] +
                            ui[1:-1,  :-2] *self.abovesea3[1:-1, 1:-1] +
                            ui[1:-1, 2:]   *self.abovesea4[1:-1, 1:-1]) /
                          self.neighboursbelow[1:-1, 1:-1] )
        u[self.neighboursbelowmask] = ud[self.neighboursbelowmask]
        return u


      nci   = nc.Dataset(self.infile,  'r')
      # take out last 3 longitude values, its double data in Brios
      print nci.variables[diffuse_var][:].shape
      try:
        data_in  = np.ma.masked_invalid(np.squeeze(nci.variables[diffuse_var][tstep,0,:,0:-3]))
      except ValueError:
        data_in  = np.ma.masked_invalid(np.squeeze(nci.variables[diffuse_var][tstep,:,0:-3]))
      print data_in.shape
      nci.close()
      projdata = extend_interp(data_in)

      notmask = ~projdata.mask
      self.notmask = notmask
      # get rid of mask, diffuse cannot handle it
      ui  = np.copy(projdata)
      # set regions that are diffused to to missval, acts as initial guess
      ui[projdata.mask] = self.diffuse_missvals[diffuse_var]
      ### set all data where ice is grounded to missval
      #ui[~self.nolandmask]  = diffuse_missvals[diffuse_var]
      u   = np.copy(ui)
      ud  = np.zeros(ui.shape)

      m=0
      while m < self.diffu_t:
        if m % 100 == 0:
          print "Diffuse for m = ", m, " for data timestep ", tstep
        u  = diffuse(ui, projdata)
        ui = u
        m += 1
      return ma.array(u, mask = projdata.mask)


    diffu_data = {"tstep":tstep}

    for diffuse_var in self.diffuse_vars:
      print diffuse_var
      diffu_data[diffuse_var] = run_diffuse(diffuse_var)

    return diffu_data


  def writeNetcdf(self, diffu_data, lite):

    ncdata = {}
    for diffu_var in self.diffuse_vars:
      ncdata[diffu_var] = ma.zeros([self.nctimesteps,len(self.x),len(self.y)])

    ## collect data
    for entry in diffu_data:
      for diffu_var in self.diffuse_vars:
        ncdata[diffu_var][entry["tstep"],:,:] = entry[diffu_var]

    nci = nc.Dataset(self.infile,  'r')
    outfile = self.outfile.strip(".nc") + "_lite.nc" if lite else self.outfile

    print "create netcdf file\n" + outfile
    ncout = nc.Dataset(self.outfile, 'w', format='NETCDF3_CLASSIC')
    ncout.createDimension('time',size=None)
    ncout.createDimension('x',size=len(self.x))
    ncout.createDimension('y',size=len(self.y))
    #ncvart  = ncout.createVariable( 'thetao','float32',('time','y','x') )
    #ncvars = ncout.createVariable( 'salinity','float32',('time','y','x')  )
    #ncvartr = ncout.createVariable( 'thetao_raw','float32',('time','y','x')  )
    #ncvarm = ncout.createVariable( 'ismelt','float32',('time','y','x')  )
    nct   = ncout.createVariable( 'time','float32',('time',) )
    ncx   = ncout.createVariable( 'x','float32',('x',) )
    ncy   = ncout.createVariable( 'y','float32',('y',) )
    nclon = ncout.createVariable( 'lon','float32',('y','x') )
    nclat = ncout.createVariable( 'lat','float32',('y','x') )

    if not lite:
      ncthk = ncout.createVariable( 'thk','float32',('y','x') )
      nctg  = ncout.createVariable( 'topg','float32',('y','x') )
      ncx1  = ncout.createVariable( 'abovesea1','float32',('y','x') )
      ncx2  = ncout.createVariable( 'abovesea2','float32',('y','x') )
      ncx3  = ncout.createVariable( 'abovesea3','float32',('y','x') )
      ncx4  = ncout.createVariable( 'abovesea4','float32',('y','x') )
      ncx5  = ncout.createVariable( 'neighboursbelowmask','float32',('y','x') )
      ncx6  = ncout.createVariable( 'neighboursbelow','float32',('y','x') )
      ncx7  = ncout.createVariable( 'nolandmask','float32',('y','x') )

    nct[:]     = self.time
    #ncvart[:]  = self.dfutemp
    #ncvars[:]  = self.dfusalt
    #ncvartr[:] = self.projtemp
    #ncvarm[:] = self.projmelt
    ncx[:]     = self.x
    ncy[:]     = self.y
    nclat[:]   = self.pismlat
    nclon[:]   = self.pismlon

    for diffu_var in self.diffuse_vars:
      ncvar    = ncout.createVariable( diffu_var,'float32',('time','y','x') )
      ncvar[:] = ncdata[diffu_var][:]
      ncvar.units = nci.variables[diffu_var].units

    if not lite:
      ncthk[:]    = self.thk
      nctg[:]     = self.topg
      ncx1[:]     = self.abovesea1
      ncx2[:]     = self.abovesea2
      ncx3[:]     = self.abovesea3
      ncx4[:]     = self.abovesea4
      ncx5[:]     = self.neighboursbelowmask
      ncx6[:]     = self.neighboursbelow
      ncx7[:]     = self.notmask
      ncthk.units  = 'meters'
      nctg.units   = 'meters'

    ncy.units = 'meters'
    ncx.units = 'meters'
    #ncvart.units = 'Kelvin'
    #ncvartr.units = 'Kelvin'
    #ncvarm.units = 'm/year'
    #ncvars.units = 'g/kg'
    nct.units    = self.timeunits
    nct.calendar = self.calendar

    ncout.diffuse_t = "diffused over " + str(self.diffu_t) + "."
    ncout.datafile = self.infile
    ncout.destination_grid_file = self.destination_grid_file
    now = datetime.datetime.now().strftime("%B %d, %Y")
    ncout.comment  = "created by matthias.mengel@pik at " + now
    ncout.history  = self.history
    nci.close()
    ncout.close()

