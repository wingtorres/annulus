import os
import sys
from inspect import getsource
from operator import attrgetter

import xarray as xr
import xgcm
import numpy as np
import math
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
#import cmocean
from datetime import timedelta as delta
import cftime

from pycoawst.tools.grid import sigma2z, metricize
from pycoawst.tools.circulation import streamwise_normal, lagrangianVelocity
from pycoawst.tools.momentum import momn2swnm

from parcels import AdvectionRK4, ErrorCode, FieldSet, JITParticle, ParticleFile, ParticleSet, Field, Variable

from matplotlib import rc
rc('text', usetex=False) 

# from dask.distributed import Client
# Client()

home = os.environ['HOME']

## Annulus
grdname = home + '/Dropbox/Oceanography/Projects/annulus/grid/ann_grid.nc'
avgname = home + '/Dropbox/Oceanography/Projects/annulus/temp/full_avg6hr_-30.nc'
dianame =  home + '/Dropbox/Oceanography/Projects/annulus/temp/full_dia6hr_-30.nc'
os.chdir(home + '/Dropbox/Oceanography/Projects/annulus/floats/')
## 1/12 Annulus
# gridname = home + '/Dropbox/Oceanography/Projects/annulus/grid/ann_grid.nc'
# filename = home + '/Desktop/anntest_strided2D.nc'
# os.chdir(home + '/Dropbox/Oceanography/Projects/annulus/floats/')

dg = xr.open_mfdataset(grdname, decode_times = False, combine = "by_coords", parallel = True)
dc = xr.open_mfdataset(avgname, decode_times = False, combine = "by_coords", parallel = True)
dm = xr.open_mfdataset(dianame, decode_times = False, combine = "by_coords", parallel = True)
dc = dc.rename({'eta_u': 'eta_rho', 'xi_v': 'xi_rho', 'xi_psi': 'xi_u', 'eta_psi': 'eta_v'})
dm = dm.rename({'eta_u': 'eta_rho', 'xi_v': 'xi_rho', 'xi_psi': 'xi_u', 'eta_psi': 'eta_v'})
sigma2z(dc, vcoord = 's_rho')
sigma2z(dc, vcoord = 's_w')

coords={'xi':{'center':'xi_rho', 'inner':'xi_u'}, 
        'eta':{'center':'eta_rho', 'inner':'eta_v'}, 
        's':{'center':'s_rho', 'outer':'s_w'}}

grid = xgcm.Grid(dc, coords = coords, periodic = 'xi')
metric = metricize(dc, grid = grid)
grid = xgcm.Grid(dc, coords = coords, metrics = metric, periodic = 'xi') #create 3D grid w/ metrics
streamwise_normal(dc, grid = grid) #rotate streamwise normal velocity field
lagrangianVelocity(dc, grid = grid) #fetch lagrangian velocity
__, dsw, dnm = momn2swnm(dm, grid) #convert to streamwise normal momentum budget

############################

# lon = Ds['lon_psi'].values
# lat = Ds['lat_psi'].values
x = dg['x_psi'].values
y = dg['y_psi'].values
time = dc['ocean_time']

#transform grid to polar
r = 6378e3 #radius of earth
z = (r**2 - x**2 - y**2)**0.5
lon, lat = (180/np.pi)*np.arctan2(y,x) , 90-(180/np.pi)*np.arccos( z/r )

#Extract lagrangian velocity from dataset
u = -dc.u_lagrangian_psi.isel(s_rho = -1).values
v = -dc.v_lagrangian_psi.isel(s_rho = -1).values

#Create parcels field set
data = {'U': u, 'V': v} 
dimensions = {'lon': lon, 'lat': lat, 'time': time}
fieldset = FieldSet.from_data(data, dimensions, transpose = False, mesh = 'spherical', allow_time_extrapolation = True) #, time_periodic = time[-1] - time[0])
fieldset.add_constant('halo_west', fieldset.U.grid.lon[0])
fieldset.add_constant('halo_east', fieldset.U.grid.lon[-1])
fieldset.add_periodic_halo(zonal = True)

h = grid.interp(grid.interp(dc.h, 'eta'),'xi')
depthField = Field(name = "h", data = h.values, lon = lon, lat = lat, time = time.values, transpose = False, mesh = 'spherical', allow_time_extrapolation = True)
fieldset.add_field(depthField)

varlist = ["accel_bar", "cor_bar", "fsco_bar", "hadv_bar", 
           "hjvf_bar", "hvisc_bar", "kvrf_bar", "prsgrd_bar",
           "bstr_bar", "wbrk_bar"]

#streamwise/normal choice here
for var in varlist:
    varval = dsw[var].values #.isel(s_rho = -1).values #3D
    field = Field(name = var, data = varval, lon = lon, lat = lat, time = time.values, transpose = False, mesh = 'spherical', allow_time_extrapolation = True)
    fieldset.add_field(field)

varlist.append("h")

#new particle class to fetch momentum budget info
class momnParticle(JITParticle):
    for v in varlist:
        setattr(JITParticle, v, Variable(v, dtype = np.float32, initial = attrgetter(v)) )

def Sample(particle, fieldset, time):
    particle.h = fieldset.h[time, particle.depth, particle.lat, particle.lon]
    particle.accel_bar = fieldset.accel_bar[time, particle.depth, particle.lat, particle.lon]
    particle.cor_bar = fieldset.cor_bar[time, particle.depth, particle.lat, particle.lon]
    particle.fsco_bar = fieldset.fsco_bar[time, particle.depth, particle.lat, particle.lon]
    particle.hadv_bar = fieldset.hadv_bar[time, particle.depth, particle.lat, particle.lon]
    particle.hjvf_bar = fieldset.hjvf_bar[time, particle.depth, particle.lat, particle.lon]
    particle.hvisc_bar = fieldset.hvisc_bar[time, particle.depth, particle.lat, particle.lon]
    particle.kvrf_bar = fieldset.kvrf_bar[time, particle.depth, particle.lat, particle.lon]
    particle.prsgrd_bar = fieldset.prsgrd_bar[time, particle.depth, particle.lat, particle.lon]
    particle.bstr_bar = fieldset.bstr_bar[time, particle.depth, particle.lat, particle.lon]
    particle.wbrk_bar = fieldset.wbrk_bar[time, particle.depth, particle.lat, particle.lon]

# def Sample(particle, fieldset, time):
#     pass

# for v in varlist:
#     setattr(Sample, v, Variable(v, dtype = np.float32, initial = attrgetter(v)) )

# print(dnm.hadv_bar.values.shape)
# plt.pcolormesh(lon,lat, np.squeeze(dnm.hadv_bar.values), vmin = -1e-5, vmax = 1e-5, cmap = "bwr")
# plt.colorbar()
# plt.show()
# assert False

#recovery kernel
def DeleteParticle(particle, fieldset, time):
    print("Deleting particle")
    particle.delete()

#periodic BC kernel
def periodicBC(particle, fieldset, time):
    if particle.lon < fieldset.halo_west:
        particle.lon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle.lon -= fieldset.halo_east - fieldset.halo_west

recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
            ErrorCode.ErrorThroughSurface: DeleteParticle}

nt = 2**4
nr = 2**2
lonp, latp = np.meshgrid( np.linspace(-180, 180, nt), np.linspace( 89.885, 89.909, nr ) )

#lonp, latp  = np.meshgrid( np.linspace(13,17,10), 89.891)
pset = ParticleSet.from_list(fieldset = fieldset, pclass = momnParticle, time = time[0], lon = lonp.flatten().tolist(), lat = latp.flatten().tolist() )

kernels = AdvectionRK4 + pset.Kernel(Sample) #+ pset.Kernel(periodicBC) 
output_file = pset.ParticleFile(name= "annFloats", outputdt = delta(seconds = 300) )
pset.execute( kernels, runtime = delta(hours = 48.0), dt = delta(seconds = 60), output_file = output_file, recovery = recovery ) #48 hours for paper
output_file.export()
output_file.close()

#sandbox
# lon, lat = .1 * 
#np.meshgrid( np.linspace(500,2500,n), np.linspace(50,2050,n) )
#pset = ParticleSet(fieldset = fieldset, lon = lon.flatten(), lat = lat.flatten(), depth = [-.25]*n*n, pclass = JITParticle)
#pset = ParticleSet.from_line(fieldset = fieldset, size = 10, start = (0, .1), finish = (0,.15), pclass = JITParticle)
# class customParticle(JITParticle):
#     theta = Variable('theta', dtype=np.float32, initial=0.)
#     r = Variable('r', dtype=np.float32, initial=0.)
    #uveitheta = Variable('uveitheta', dtype = np.complex128, initial=0.)

#    # particle.theta = math.atan2(particle.lat, particle.lon)
    # particle.r = math.hypot(particle.lon, particle.lat)

    # if particle.theta >= math.pi/12:
    #     # particle.uveitheta = particle.r*math.exp(1j*(particle.theta - math.pi/6 ))
    #     particle.lon = particle.r*math.cos(particle.theta - math.pi/6) #uveitheta.real
    #     particle.lat = particle.r*math.sin(particle.theta - math.pi/6) #uveitheta.imag

    # elif particle.theta <= -math.pi/12:
    #     # particle.uveitheta = particle.r*math.exp(1j*(particle.theta + math.pi/6 ))
    #     particle.lon = particle.r*math.cos(particle.theta + math.pi/6) #uveitheta.real
    #     particle.lat = particle.r*math.sin(particle.theta + math.pi/6) #uveitheta.imag
    #sanity check on velocity field
# plt.figure()
# plt.quiver( Ds.lon_psi, Ds.lat_psi, u_eastward.isel(ocean_time = 0), v_northward.isel(ocean_time = 0) )
# plt.pcolormesh(Ds.lon_psi, Ds.lat_psi, v_northward.isel(ocean_time = 0), cmap = cmocean.cm.balance, vmin = -0.5, vmax = 0.5)
# plt.show()
# sys.exit()
# pset = ParticleSet.from_line(fieldset = fieldset, size = 10, start = (.095, 0), finish = (.15,0), pclass = JITParticle )
# lonp = .1*np.cos(theta)
# latp = .1*np.sin(theta)