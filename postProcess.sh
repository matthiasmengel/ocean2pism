#!/bin/bash

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

destdir=`echo $PWD | awk -F/ '{print $NF}'`
#sourcepath=/iplex/01/tumble/mengel/pismSourceData/20121113_RacmoEcEarthTemperatureAndSmbCorrected
#destpath=/iplex/01/tumble/mengel/pismInputData/$destdir/multiyearave
destpath=/scratch/01/mengel/pismInputData/$destdir
mkdir -p $destpath/multiyearave

#create 30 year time mean of historical data
ncra -F -O -d time,1,30,1 $destpath/diffused/BRIOS.HadGem2_20C_sigma24_diffused10000_lite.nc $destpath/multiyearave/BRIOS.HadGem2_20C_sigma24_diffused10000_lite_1950_1980.nc
#get rid of time dimension
ncwa -F -O -a time $destpath/multiyearave/BRIOS.HadGem2_20C_sigma24_diffused10000_lite_1950_1980.nc $destpath/multiyearave/BRIOS.HadGem2_20C_sigma24_diffused10000_lite_1950_1980.nc

# merge the hist and rcp datasets
#ncrcat -F -O $destpath/pism_15km_hist.nc $destpath/pism_15km_RCP45.nc $destpath/merged/pism_15km_hist_RCP45.nc
#ncrcat -F -O $destpath/pism_15km_hist.nc $destpath/pism_15km_RCP85.nc $destpath/merged/pism_15km_hist_RCP85.nc

 ##create anomalies of merged runs
#for scen in RCP45 RCP85 ; do

  ##ncwa -F -O -a time -d time,1 $destpath/pism_15km_$scen.nc $destpath/pism_15km_${scen}_1styear.nc
  #ncdiff -O $destpath/merged/pism_15km_hist_$scen.nc $destpath/racmoEC_EARTH_15km_hist_19501980.nc $destpath/merged/pism_15km_hist_${scen}_anom.nc
  #ncrename -O -v ice_surface_temp,ice_surface_temp_anomaly -v climatic_mass_balance,climatic_mass_balance_anomaly $destpath/merged/pism_15km_hist_${scen}_anom.nc $destpath/merged/pism_15km_hist_${scen}_anom.nc
  #cpy rewriteTimeBounds.py $destpath/merged/pism_15km_hist_${scen}_anom.nc
  #ncap2 -O -s 'climatic_mass_balance_anomaly=climatic_mass_balance_anomaly*0' $destpath/merged/pism_15km_hist_${scen}_anom.nc $destpath/merged/pism_15km_hist_${scen}_anom_nosmb.nc
  #ncap2 -O -s 'ice_surface_temp_anomaly=ice_surface_temp_anomaly*0' $destpath/merged/pism_15km_hist_${scen}_anom.nc $destpath/merged/pism_15km_hist_${scen}_anom_notemp.nc

#done
