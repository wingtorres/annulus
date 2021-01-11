#!/usr/bin/env python
# coding: utf-8

import os
import sys
import panel as pn
import param
import numpy as np
from scipy.interpolate import RegularGridInterpolator as rgi

import holoviews as hv
import hvplot.xarray
import hvplot.pandas
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import rc
import cmocean

from parcels import AdvectionRK4, AdvectionRK45, ErrorCode, Variable, Field, FieldSet, JITParticle, ParticleFile, ParticleSet
from datetime import timedelta as delta

import pycoawst.tools.grid as pcg
import pycoawst.tools.circulation as pcc
import pycoawst.tools.momentum as pcm

pn.extension()
hv.extension('bokeh')
home = os.environ["HOME"]
rc("text", usetex=False)

from dask.distributed import Client
Client()

#display(sys.argv[2])
#avgname =  f"{sys.argv[1]}/results/ocean_avg_00239.nc"
#dianame =  f"{sys.argv[1]}/results/ocean_dia_00239.nc"
#cd = sys.argv[2]

ds = xr.open_dataset(avgname, use_cftime = True).isel(ocean_time = -1)
ds, grid = pcg.xromsgrid(ds, vertical = False)
ds = pcc.strain_tensor(ds, grid)
ds = pcc.streamwise_normal(ds, grid, vertical = False)
ds = pcc.velocity_rho(ds,grid)
ds = pcc.shear_orbital(ds)
#ds = pcc.curvature_sn(ds, grid, .025)
ds = pcc.volume_budget(ds, grid)
#ds = pcc.pressure_torque(ds, grid)

dm = xr.open_dataset(dianame, decode_times = False, use_cftime = True).isel(ocean_time = -1)
dm = pcg.recoord(dm)
dm["pm_psi"], dm["pn_psi"] = ds.pm_psi, ds.pn_psi
dm, dsw, dnm = pcm.momn2swnm(dm, grid)

dm = pcc.strain_tensor_psi(dm, grid)
dm = pcc.shear_orbital(dm)
dm = pcc.vorticity(dm, ds, grid, dsw, dnm)
dm = pcc.curvature_dia(dm, ds, dnm, grid, cd = cd)

class witParticle(JITParticle):
    dvdn = Variable('dvdn',  dtype=np.float32, initial=0.)
    dvdn_int = Variable('dvdn_int', dtype=np.float32, initial=0.)

def Sample(particle, fieldset, time):
    particle.dvdn = fieldset.dvdn[time, particle.depth, particle.lat, particle.lon]
    ds = (particle.lon**2+ particle.lat**2)**0.5
    particle.dvdn_int += abs(particle.dvdn)*ds #integrate absolute value

def DeleteParticle(particle, fieldset, time):
    print("Deleting particle")
    particle.delete()

recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
            ErrorCode.ErrorThroughSurface: DeleteParticle}

ds["lon_psi"], ds["lat_psi"] = pcg.cart2polar(ds.y_psi, ds.x_psi, pole="south", phase = -90)
ds = ds.set_coords(["lon_psi","lat_psi"])

dm["lon_psi"], dm["lat_psi"] = pcg.cart2polar(dm.y_psi, dm.x_psi, pole="south", phase = -90)
dm = dm.set_coords(["lon_psi","lat_psi"])

dm["anglep"] = np.arctan2(dm.y_psi, dm.x_psi) - np.pi/2
dm["ubar_eastward"], dm["vbar_northward"] = pcc.curv2cart(dm.ubar_psi, dm.vbar_psi, dm.anglep)

velocities = {'U': 'ubar_eastward',
             'V': 'vbar_northward'}

dimensions = {'U': {'lon': 'x_psi', 'lat': 'y_psi'},
              'V': {'lon': 'x_psi', 'lat': 'y_psi'}}

mesh = "flat"
x = dm.x_psi
y = dm.y_psi
fieldset = FieldSet.from_xarray_dataset(dm, variables = velocities, dimensions = dimensions, mesh = mesh, time_periodic = False)
field_dvdn = Field(name = "dvdn", data = dm.dvdn.values, lon = x, lat = y, transpose = False, mesh = mesh, allow_time_extrapolation = True)
fieldset.add_field(field_dvdn)

def find_trajectory(phi, ds, r = 12.05e3):
    """Find trajectory that minimizes path integral of the given variable"""

    x0,y0 = r*np.cos(phi), r*np.sin(phi)

    pset = ParticleSet.from_list(fieldset = fieldset, pclass = witParticle, time = ds.ocean_time.values, lon = x0, lat = y0 )
    #output_file = pset.ParticleFile(name="output.nc", outputdt=delta(seconds=60))
    kernels = AdvectionRK45 + pset.Kernel(Sample)
    pset.execute(kernels, runtime = delta(hours = 6.0), dt = delta(seconds = 60), recovery = recovery)#, output_file=output_file)
    return abs(pset.dvdn_int[0])


from scipy.optimize import minimize

options={"disp": True}
F = minimize(find_trajectory, 0, args=dm , method = None, tol=1e-3, options = options, bounds = [ (-np.pi/192, np.pi/192) ] )
F

#np.save(f"{sys.argv[1]}/theta_opt", F.x)
np.save(f"{resultsdir}/theta_opt", F.x)