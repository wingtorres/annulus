import os
import sys
import math
import numpy as np
import xarray as xr

from datetime import timedelta as delta
import cftime

from parcels import AdvectionRK4, ErrorCode, FieldSet, JITParticle, ParticleFile, ParticleSet, Variable, plotting, Kernel

import romspy
import romspy.accessor

home = os.environ['HOME']

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

def roms2tracks(gridname, filename, outfile = "annfloats_slice.nc"):
    
    print(gridname, filename, outfile)
    ds = xr.open_dataset(filename, decode_times = False)
    dg = xr.open_dataset(gridname)

    time = ds['ocean_time']
    x = dg['x_psi'].values
    y = dg['y_psi'].values
    r,t = np.arctan2(y,x), np.hypot(x,y)

    ds.xwit.make_grid
    ds.xwit.lagrangian_velocity(rotate = True)

    #data = {'U': -ds.xwit.ds.ubar_lagrangian_psi.values, 'V': -ds.xwit.ds.vbar_lagrangian_psi.values} #2D velocity
    data = {'U': ds.xwit.ds.vbar_lagrangian_psi.values, 'V': -ds.xwit.ds.ubar_lagrangian_psi.values} #2D velocity cartesian
    dimensions = {'lon': x, 'lat': y}
    fieldset = FieldSet.from_data(data, dimensions, transpose = False, mesh = "flat")

    recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
                ErrorCode.ErrorThroughSurface: DeleteParticle}

    nt = 2**6
    nr = 2**5
    rp, tp = np.meshgrid( np.linspace(10.5e3, 13.5e3, nr), np.linspace(-14*np.pi/180 , 14*np.pi/180, nt) ) 
    xp, yp = rp*np.cos(tp), rp*np.sin(tp)
    xp, yp = xp.flatten().tolist(), yp.flatten().tolist()

    pset = ParticleSet.from_list(fieldset = fieldset, pclass = witParticle, lon = xp, lat = yp )

    # fieldset.U.show()
    # plt.pcolormesh(x, y, dg.h.values[:-1,:-1])
    # plt.scatter(xp,yp, s = 10, color = "k")
    # plt.show()
    #assert False

    kernels = AdvectionRK4 + pset.Kernel(PeriodicParticle)
    output_file = pset.ParticleFile(name = outfile, outputdt = delta(seconds = 30) )
    pset.execute( kernels, runtime = delta(hours = 60.0), dt = delta(seconds = 30), output_file = output_file, recovery = recovery)
    output_file.export()
    output_file.close()


def main():
    gridname = sys.argv[1]
    avgname = sys.argv[2]
    outname = sys.argv[3]
    roms2tracks(gridname, avgname, outname)

if __name__ == "__main__":
    main()
