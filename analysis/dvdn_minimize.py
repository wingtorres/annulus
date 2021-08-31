#!/usr/bin/env python
# coding: utf-8

import os
import sys
import math
import numpy as np
from scipy.optimize import minimize, minimize_scalar, differential_evolution
import xarray as xr

import matplotlib.pyplot as plt
import cmocean

from parcels import AdvectionRK4, AdvectionRK45, ErrorCode, Variable, Field, FieldSet, JITParticle, ParticleFile, ParticleSet
from datetime import timedelta as delta
from operator import attrgetter

# import pycoawst.tools.grid as pcg
# import pycoawst.tools.circulation as pcc
# import pycoawst.tools.momentum as pcm

# home = os.environ["HOME"]

#todo:
#-stop integration after certain distance traveled

class witParticle(JITParticle):
    s = Variable('s', dtype = np.float32, initial=0.)
    lon_lag = Variable("lon_lag", dtype = np.float32, to_write = False, initial = attrgetter("lon"))
    lat_lag = Variable("lat_lag", dtype = np.float32, to_write = False, initial = attrgetter("lat"))
    #var_lag = Variable("lon_lag", dtype = np.float32, to_write = False, initial = attrgetter(dvdn))
    #dvdn = Variable('dvdn',  dtype=np.float32, initial=attrgetter("dvdn"))
    #dvdn_int = Variable('dvdn_int', dtype=np.float32, initial=0.)
    #dadn = Variable('dadn',  dtype=np.float32, initial=attrgetter("dadn"))
    #dadn_int = Variable('dadn_int',  dtype=np.float32, initial=0.)

    shearobj = Variable('shearobj', dtype=np.float32, initial=attrgetter("shearobj"))
    shearobj_int = Variable('shearobj_int', dtype=np.float32, initial=0.)


def objective_function(dvdn, v = 1, dadn = 0):
    return abs(dvdn/v) + abs(dadn)

def Sample(particle, fieldset, time):
    # if particle.s == particle.s:
    dx = particle.lon - particle.lon_lag
    dy = particle.lat - particle.lat_lag
    ds = math.sqrt( math.pow(dx, 2) + math.pow(dy, 2) )
    particle.s += ds
    # particle.dvdn = fieldset.dvdn[time, particle.depth, particle.lat, particle.lon]
    # particle.dvdn_int += math.fabs(particle.dvdn)*ds #integrate absolute value
    #particle.dadn = fieldset.dadn[time, particle.depth, particle.lat, particle.lon]
    #particle.dadn_int += math.fabs(particle.dadn)*ds
    
    particle.shearobj = math.fabs( fieldset.dvdn[time, particle.depth, particle.lat, particle.lon] )#/fieldset.v[time, particle.depth, particle.lat, particle.lon] )
    particle.shearobj += math.fabs( fieldset.v[time, particle.depth, particle.lat, particle.lon]*fieldset.dadn[time, particle.depth, particle.lat, particle.lon] )
    
    #particle.shearobj = math.fabs( fieldset.dvdn[time, particle.depth, particle.lat, particle.lon] )
    #particle.shearobj = abs(fieldset.dvdn[time, particle.depth, particle.lat, particle.lon]/fieldset.v[time, particle.depth, particle.lat, particle.lon])
    #particle.shearobj = math.fabs(fieldset.dadn[time, particle.depth, particle.lat, particle.lon])

    particle.shearobj_int += particle.shearobj*ds #/max(particle.s,1) #integrate absolute value, weight by distance along track
    
    particle.lon_lag = particle.lon
    particle.lat_lag = particle.lat

    if particle.s >= 4e3:
        particle.state = 4
        print("Halting execution: Max distance reached")

def DeleteParticle(particle, fieldset, time):
    particle.state = 4 
    print("Halting Execution: Out of bounds")
    #ErrorCode.Success
    #print(f"Deleting particle...s={particle.s:0.2f} m")
    #particle.delete()

def find_trajectory(x, ds, pset, t0, scale, r = 12.0125e3):
    """Find trajectory that minimizes path integral of the given variable"""
    #phi = x[0] #for general minimization algorithm
    phi = x/scale #for minimize_scalar (scaled)
    x0, y0 = r*np.cos(phi), r*np.sin(phi)

    #display(phi, x0, y0)
    recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
                ErrorCode.ErrorThroughSurface: DeleteParticle}
    
    
#     pset._collection.data["time"] = t0
#     pset._collection.data["shearobj_int"] = np.array( [0.], dtype="float32")
#     pset._collection.data["s"] = np.array( [0.], dtype="float32")
#     pset._collection.data["lon"], pset._collection.data["lat"] = np.array( [x0], dtype="float32"),  np.array( [y0], dtype="float32")
#     pset._collection.data["lon_lag"], pset._collection.data["lat_lag"] = np.array( [x0], dtype="float32"),  np.array( [y0], dtype="float32")
    
#    display(pset.particle_data)
    pset.particle_data["time"] = t0
    pset.particle_data["shearobj_int"] = np.array( [0.], dtype="float32")
    pset.particle_data["s"] = np.array( [0.], dtype="float32")
    pset.particle_data["lon"], pset.particle_data["lat"] = np.array( [x0], dtype="float32"),  np.array( [y0], dtype="float32")
    pset.particle_data["lon_lag"], pset.particle_data["lat_lag"] = np.array( [x0], dtype="float32"),  np.array( [y0], dtype="float32")

    kernels = AdvectionRK4 + pset.Kernel(Sample) #+ pset.Kernel(Stop)
    pset.execute(kernels, runtime = delta(hours = 120.0), dt = delta(seconds = 30), recovery = recovery, verbose_progress = True)
    objective = pset.shearobj_int[0]/pset.s[0] #min(pset.s[0], 1500) #pset.s[0]
    #objective =  math.fabs(pset.dvdn_int[0]) # + math.fabs(pset.vdadn_int[0])
    
    return objective

def optimize_trajectory(dmpath, outfile):
    dm = xr.load_dataset(dmpath, decode_times = False, use_cftime = True)

    velocities = {'U': 'ubar_eastward',
                  'V': 'vbar_northward'}

    dimensions = {'U': {'lon': 'x_psi', 'lat': 'y_psi'},
                  'V': {'lon': 'x_psi', 'lat': 'y_psi'}}

    mesh = "flat"
    x,y = dm.x_psi, dm.y_psi
    fieldset = FieldSet.from_xarray_dataset(dm, variables = velocities, dimensions = dimensions, mesh = mesh, time_periodic = False)
    field_dvdn = Field(name = "dvdn", data = dm.dvdn.values, lon = x, lat = y, transpose = False, mesh = mesh, allow_time_extrapolation = True)
    field_dadn = Field(name = "dadn", data = dm.dadn.values, lon = x, lat = y, transpose = False, mesh = mesh, allow_time_extrapolation = True)
    field_v = Field(name = "v", data = dm.V.values, lon = x, lat = y, transpose = False, mesh = mesh, allow_time_extrapolation = True)
    fieldset.add_field(field_dvdn)
    fieldset.add_field(field_dadn)
    fieldset.add_field(field_v)
    #p = dm.dadn.plot(x = "x_psi", y = "y_psi", cmap = "twilight_shifted", vmin = -1e-3, vmax = 1e-3)
    #plt.show(p)

    pset = ParticleSet.from_list(fieldset = fieldset, pclass = witParticle, time = dm.ocean_time.values, lon = 0, lat = 0)
    #display(pset["time"])

    
    t0 = pset.particle_data["time"] #._collection.data["time"]

    scale = 1e3
    options={"disp": True, "xatol": 1e-7}
    bounds = ( scale*(-np.pi/384 + np.pi/2), scale*(np.pi/384 + np.pi/2) )
    F = minimize_scalar(find_trajectory, scale*np.pi/2, args=(dm,pset,t0,scale), method = "bounded", options = options, bounds = bounds)
    #F = minimize(find_trajectory, 0, args=dm, method = "SLSQP", tol = None, options = options, bounds = [ (-np.pi/192, np.pi/192) ] )
    #F = differential_evolution(find_trajectory, args=(dm,pset,t0,scale,),  disp = True, bounds = [bounds])# , workers = -1)
    
    display(F)
    np.save(outfile, F.x/scale)
    dm.close()

