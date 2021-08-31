import os
from operator import attrgetter
import sys
import math
import numpy as np
import xarray as xr

from datetime import timedelta as delta
import cftime

import matplotlib.pyplot as plt

from parcels import AdvectionRK4, ErrorCode, FieldSet, JITParticle, ParticleFile, ParticleSet, Variable, plotting, Kernel

import romspy
import romspy.accessor
from romspy.grid import recoord

home = os.environ['HOME']

class witParticle(JITParticle):
    theta = Variable('theta', dtype=np.float32, initial=0.)
    r = Variable('r', dtype=np.float32, initial=0.)
    s = Variable('s', dtype = np.float32, initial=0.)
    lon_lag = Variable("lon_lag", dtype = np.float32, to_write = False, initial = attrgetter("lon"))
    lat_lag = Variable("lat_lag", dtype = np.float32, to_write = False, initial = attrgetter("lat"))
        
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

def PathDistance(particle, fieldset, time):
    dx = particle.lon - particle.lon_lag
    dy = particle.lat - particle.lat_lag
    ds = math.sqrt( math.pow(dx, 2) + math.pow(dy, 2) )
    particle.s += ds
    
    particle.lon_lag = particle.lon
    particle.lat_lag = particle.lat
    
    if particle.s >= 4e3:
        particle.state = 4
        particle.delete()
        #print("Halting execution: Max distance reached")
    
def DeleteParticle(particle, fieldset, time):
    particle.state = 4
    print("Deleting particle")
    particle.delete()

def strain_ev(ds, parameter = 1.0):
    a = ds.dudx
    b = 0.5*(ds.dvdx + ds.dudy)
    c = 0.5*(ds.dvdx + ds.dudy)
    d = ds.dvdy
    
    ds["trace"] = ds.dudx + ds.dvdy
    ds["det"] = ds.dudx*ds.dvdy - 0.25*(ds.dvdx + ds.dudy)**2
    ds["ev1"] = 0.5*(ds.trace + np.sqrt(ds.trace**2 - 4*ds.det))
    ds["ev2"] = 0.5*(ds.trace - np.sqrt(ds.trace**2 - 4*ds.det))
    
    #http://math.colgate.edu/~wweckesser/math312Spring06/handouts/IMM_2x2linalg.pdf
    ds["ev1_x"] = b 
    ds["ev1_y"] = ds.ev1 - a
    ds["ev2_x"] = ds.ev2 - d
    ds["ev2_y"] = c
    
    ds["eta1_x"] = np.sqrt( (+ds.ev2 - parameter**2)/(ds.ev2 - ds.ev1) + 0j)*ds.ev1_x \
                 + np.sqrt( (-ds.ev1 + parameter**2)/(ds.ev2 - ds.ev1) + 0j)*ds.ev2_x
    ds["eta1_y"] = np.sqrt( (+ds.ev2 - parameter**2)/(ds.ev2 - ds.ev1) + 0j)*ds.ev1_y \
                 + np.sqrt( (-ds.ev1 + parameter**2)/(ds.ev2 - ds.ev1) + 0j)*ds.ev2_y
    ds["eta2_x"] = np.sqrt( (+ds.ev2 - parameter**2)/(ds.ev2 - ds.ev1) + 0j)*ds.ev1_x \
                 - np.sqrt( (-ds.ev1 + parameter**2)/(ds.ev2 - ds.ev1) + 0j)*ds.ev2_x
    ds["eta2_y"] = np.sqrt( (+ds.ev2 - parameter**2)/(ds.ev2 - ds.ev1) + 0j)*ds.ev1_y \
                 - np.sqrt( (-ds.ev1 + parameter**2)/(ds.ev2 - ds.ev1) + 0j)*ds.ev2_y
    
    ds["eta1_x"] = np.real(ds.eta1_x)
    ds["eta1_y"] = np.real(ds.eta1_y)
    ds["eta2_x"] = np.real(ds.eta2_x)
    ds["eta2_y"] = np.real(ds.eta2_y)
    return ds
    
def roms2tracks(gridname, filename, outfile = "lcsfloats.nc", lcs = "repelling", direction = "forward"):
    
    print(gridname, filename, outfile)
    ds = xr.open_dataset(filename, decode_times = False)
    dg = xr.open_dataset(gridname)
    dg = recoord(dg)
    
    ds.xwit.make_grid
    ds.xwit.recoord
    ds.xwit.ds = ds.xwit.ds.assign_coords( {"x_psi": dg.coords["x_psi"], "y_psi": dg.coords["y_psi"]})

    ds.xwit.streamwise_normal
    #ds.xwit.lagrangian_velocity(rotate = True)
    ds.xwit.strain_tensor(uname = "ubar_eastward", vname = "vbar_northward", interp2psi = False)
    ds.xwit.ds = strain_ev(ds.xwit.ds)
        
    time = ds['ocean_time']
    x = ds.xwit.ds.coords['x_psi'].values
    y = ds.xwit.ds.coords['y_psi'].values
    r,t = np.hypot(x,y), np.arctan2(y,x)
    t = t - np.pi/2
    
    amp = .001
    if lcs == "attracting":
        u,v = ds.xwit.ds.ev1_x.values, ds.xwit.ds.ev1_y.values
    elif lcs == "repelling":
        u,v = ds.xwit.ds.ev2_x.values, ds.xwit.ds.ev2_y.values
    elif lcs == "elliptic":
        u,v = ds.xwit.ds.eta1_x.values, ds.xwit.ds.eta1_y.values
        amp = 1
        
    norm = amp*max(np.hypot(u,v).flatten())
    
    data = {'U': u/norm, 'V': v/norm} #2D velocity cartesian
    dimensions = {'lon': x, 'lat': y}
    fieldset = FieldSet.from_data(data, dimensions, transpose = False, mesh = "flat")

    recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
                ErrorCode.ErrorThroughSurface: DeleteParticle}
    
    #cn local maxima applied every cnth point
    #https://stackoverflow.com/questions/63093258/find-values-within-a-radius-from-multiples-lat-lon-centers-in-a-netcdf
    cn = 10
    local_radius = 250 #meters
    centers = np.array( list(zip( x[::cn,::cn].flatten(), y[::cn,::cn].flatten() )) )
    ds_c = xr.DataArray(centers, coords={"coord": ["x_psi", "y_psi"]}, dims = ["center", "coord"])

    ds.xwit.ds["distance"] = np.hypot( ds_c.sel(coord = "x_psi") - ds.xwit.ds.coords["x_psi"], 
                                       ds_c.sel(coord = "y_psi") - ds.xwit.ds.coords["y_psi"])
       
    if lcs == "attracting":
        ev  = abs( ds.xwit.ds.ev2.where(ds.xwit.ds.distance < local_radius) )
    elif lcs == "repelling":
        ev  = abs( ds.xwit.ds.ev1.where(ds.xwit.ds.distance < local_radius) )
    
    if lcs == "attracting" or lcs == "repelling":
        inds = ev.argmax(dim = ["eta_v", "xi_u"], skipna = True)
        xp, yp = ds.xwit.ds.coords["x_psi"].isel(inds).values, ds.xwit.ds.coords["y_psi"].isel(inds).values
    
    #check if points are within local radius of eachother then take max?
    
    if lcs == "elliptic":
        #poincare sections:
        xp1, yp1 = np.linspace(10.5e3, 14e3), np.linspace(1e3, 2e3)
        xp2, yp2 = np.linspace(10.5e3, 14e3), np.linspace(-1e3, -2e3)
        xp, yp = np.concatenate( [xp1,xp2] ), np.concatenate( [yp1,yp2] )
        
    #xp, yp = list(set(list(xp))), list(set(list(yp))) #make unique
    pset = ParticleSet.from_list(fieldset = fieldset, pclass = witParticle, lon = xp, lat = yp )

    #fieldset.U.show()
    #plt.pcolormesh(x, y, ds.xwit.ds.vbar_northward, cmap = "bwr", vmin = -.1, vmax = .1)
#     plt.pcolormesh(x, y, ds.xwit.ds.vbar_northward, cmap = "bwr", vmin = -.1, vmax = .1)
#     plt.colorbar()
#     plt.scatter(xp, yp, s = 10, color = "k")
#     plt.show()
#     assert False

    if direction == "forward":
        sign = 1
    elif direction == "backward":
        sign= -1
        
    kernels = AdvectionRK4 + pset.Kernel(PeriodicParticle) #+ pset.Kernel(PathDistance)
    output_file = pset.ParticleFile(name = f"lcsfloats_{lcs}_{direction}.nc", outputdt = delta(seconds = 30) )
    pset.execute( kernels, runtime = delta(hours = 12.0), dt = sign*delta(seconds = 30), output_file = output_file, recovery = recovery)
    output_file.export()
    output_file.close()

def main():
    gridname = sys.argv[1]
    avgname = sys.argv[2]
    outname = sys.argv[3]
    lcs = sys.argv[4]
    direction = sys.argv[5]
    roms2tracks(gridname, avgname, outname, lcs, direction)

if __name__ == "__main__":
    main()
