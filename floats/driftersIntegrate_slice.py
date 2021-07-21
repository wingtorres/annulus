import os
import sys
import xarray as xr
# import xgcm
import numpy as np
import math
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from datetime import timedelta as delta
import cftime

from parcels import AdvectionRK4, ErrorCode, FieldSet, JITParticle, ParticleFile, ParticleSet, Variable, plotting, Kernel

home = os.environ['HOME']

import romspy
import romspy.accessor
from romspy import circulation as pcc
from romspy.grid import cart2polar


## 1/12 Annulus
gridname = home + '/Dropbox/Oceanography/Projects/annulus/grid/ann_grid_smol.nc'
avgname = '../output/annulus_sector_3D/smol_lat_-30_z0_0.25/results/ocean_avg.nc'
os.chdir(home + '/Dropbox/Oceanography/Projects/annulus/floats/')

def roms2tracks(gridname, avgname, xp, yp):

ds = xr.open_dataset(filename, decode_times= False)#.isel(ocean_time = -1)
dg = xr.open_dataset(gridname)

time = ds['ocean_time']
x = dg['x_psi'].values
y = dg['y_psi'].values
r,t = np.arctan2(y,x), np.hypot(x,y)

#Transform grid to polar
lon, lat = cart2polar(x,y)

ds.xwit.make_grid
ds.xwit.lagrangian_velocity(rotate = True)

#data = {'U': -ds.xwit.ds.ubar_lagrangian_psi.values, 'V': -ds.xwit.ds.vbar_lagrangian_psi.values} #2D velocity
data = {'U': ds.xwit.ds.vbar_lagrangian_psi.values, 'V': -ds.xwit.ds.ubar_lagrangian_psi.values} #2D velocity cartesian
dimensions = {'lon': x, 'lat': y}
fieldset = FieldSet.from_data(data, dimensions, transpose = False, mesh = "flat")

class witParticle(JITParticle):
    theta = Variable('theta', dtype=np.float32, initial=0.)
    r = Variable('r', dtype=np.float32, initial=0.)
    
def PeriodicParticle(particle, fieldset, time):
    
    particle.theta = math.atan2(particle.lat, particle.lon)
    particle.r = (particle.lon**2 + particle.lat**2)**0.5
    
    dt = 29*math.pi/180
    if particle.theta >= 14.5*math.pi/180:
        #print("rotating particle")
        # particle.uveitheta = particle.r*math.exp(1j*(particle.theta - math.pi/6 ))
        particle.lon = particle.r*math.cos(particle.theta - dt) #uveitheta.real
        particle.lat = particle.r*math.sin(particle.theta - dt) #uveitheta.imag

    elif particle.theta <= -14.5*math.pi/180:
        # particle.uveitheta = particle.r*math.exp(1j*(particle.theta + math.pi/6 ))
        particle.lon = particle.r*math.cos(particle.theta + dt) #uveitheta.real
        particle.lat = particle.r*math.sin(particle.theta + dt) #uveitheta.imag

def DeleteParticle(particle, fieldset, time):
    print("Deleting particle")
    particle.delete()

recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
            ErrorCode.ErrorThroughSurface: DeleteParticle}



#lonp, latp = np.meshgrid( np.linspace(-14, 14, nt), np.linspace( 89.885, 89.909, nr ) )


lonp = lonp.flatten().tolist()
latp = latp.flatten().tolist()

xp, yp = xp.flatten().tolist(), yp.flatten().tolist()
# xp, yp = 13000, 0

pset = ParticleSet.from_list(fieldset = fieldset, pclass = witParticle, lon = xp, lat = yp )

# fieldset.U.show()
# plt.pcolormesh(x, y, dg.h.values[:-1,:-1])
# plt.scatter(xp,yp, s = 10, color = "k")
# plt.show()
#assert False

kernels = AdvectionRK4 + pset.Kernel(PeriodicParticle)
output_file = pset.ParticleFile(name= "annfloats_slice", outputdt = delta(seconds = 30) )
pset.execute( kernels, runtime = delta(hours = 60.0), dt = delta(seconds = 30), output_file = output_file, recovery = recovery )
output_file.export()
output_file.close()

#sandbox