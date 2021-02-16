#!/usr/bin/env python
# coding: utf-8

import os
import sys
import math
import numpy as np
from scipy.optimize import minimize, minimize_scalar, differential_evolution
import xarray as xr

from parcels import AdvectionRK4, AdvectionRK45, ErrorCode, Variable, Field, FieldSet, JITParticle, ParticleFile, ParticleSet
from datetime import timedelta as delta
from operator import attrgetter

# import pycoawst.tools.grid as pcg
# import pycoawst.tools.circulation as pcc
# import pycoawst.tools.momentum as pcm

# home = os.environ["HOME"]

#todo:
#-stop integration after certain distance traveled


dm = xr.load_dataset(sys.argv[1], decode_times = False, use_cftime = True)

class witParticle(JITParticle):
    s = Variable('s', dtype = np.float32, initial=0.)
    lon_lag = Variable("lon_lag", dtype = np.float32, to_write = False, initial = attrgetter("lon"))
    lat_lag = Variable("lat_lag", dtype = np.float32, to_write = False, initial = attrgetter("lat"))
    #var_lag = Variable("lon_lag", dtype = np.float32, to_write = False, initial = attrgetter(dvdn))
    dvdn = Variable('dvdn',  dtype=np.float32, initial=attrgetter("dvdn"))
    dvdn_int = Variable('dvdn_int', dtype=np.float32, initial=0.)
    
def Sample(particle, fieldset, time):
    # if particle.s == particle.s:
    dx = particle.lon - particle.lon_lag
    dy = particle.lat - particle.lat_lag
    ds = math.sqrt( math.pow(dx, 2) + math.pow(dy, 2) )
    particle.s += ds
    particle.dvdn = fieldset.dvdn[time, particle.depth, particle.lat, particle.lon]
    particle.dvdn_int += math.fabs(particle.dvdn)*ds #integrate absolute value
    particle.lon_lag = particle.lon
    particle.lat_lag = particle.lat

    if particle.s >= 4e3:
        particle.state = 4
        print("Halting execution...")

def DeleteParticle(particle, fieldset, time):
    particle.state = 4 
    print("Halting Execution: Out of bounds")
    #ErrorCode.Success
    #print(f"Deleting particle...s={particle.s:0.2f} m")
    #particle.delete()

recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
            ErrorCode.ErrorThroughSurface: DeleteParticle}

velocities = {'U': 'ubar_eastward',
              'V': 'vbar_northward'}

dimensions = {'U': {'lon': 'x_psi', 'lat': 'y_psi'},
              'V': {'lon': 'x_psi', 'lat': 'y_psi'}}

mesh = "flat"
x,y = dm.x_psi, dm.y_psi
fieldset = FieldSet.from_xarray_dataset(dm, variables = velocities, dimensions = dimensions, mesh = mesh, time_periodic = False)
field_dvdn = Field(name = "dvdn", data = dm.dvdn.values, lon = x, lat = y, transpose = False, mesh = mesh, allow_time_extrapolation = True)
fieldset.add_field(field_dvdn)

pset = ParticleSet.from_list(fieldset = fieldset, pclass = witParticle, time = dm.ocean_time.values, lon = 0, lat = 0)
t0 = pset.particle_data["time"]

def find_trajectory(x, ds, r = 12.05e3):
    """Find trajectory that minimizes path integral of the given variable"""
    #phi = x[0]
    phi = x #for minimize_scalar
    x0, y0 = r*np.cos(phi), r*np.sin(phi)

    pset.particle_data["time"] = t0
    pset.particle_data["dvdn_int"] = np.array( [0.], dtype="float32")
    pset.particle_data["s"] = np.array( [0.], dtype="float32")
    pset.particle_data["dvdn"] = np.array( [0.], dtype="float32")
    pset.particle_data["lon"], pset.particle_data["lat"] = np.array( [x0], dtype="float32"),  np.array( [y0], dtype="float32") ##np.float32(x0), np.float32(y0) 
    pset.particle_data["lon_lag"], pset.particle_data["lat_lag"] = np.array( [x0], dtype="float32"),  np.array( [y0], dtype="float32")

    kernels = AdvectionRK4 + pset.Kernel(Sample) #+ pset.Kernel(Stop)
    pset.execute(kernels, runtime = delta(hours = 1200.0), dt = delta(seconds = 60), recovery = recovery, verbose_progress = False)#, output_file=output_file)
    #print((pset))
    objective =  math.fabs(pset.dvdn_int[0])
    return objective

options={"disp": True}#, "ftol": 1e-6}#, "eps": np.pi/19200}
F = minimize_scalar(find_trajectory, 0, args=dm, method = "bounded", tol = None, options = options, bounds = (-np.pi/192, np.pi/192) )
#F = minimize(find_trajectory, 0, args=dm, method = "SLSQP", tol = None, options = options, bounds = [ (-np.pi/192, np.pi/192) ] )
#F = differential_evolution(find_trajectory, args=(dm,12.05e3),  disp = True, bounds = [ (-np.pi/192, np.pi/192) , workers = -1)

print(F)

