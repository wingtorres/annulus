import os
import numpy as np

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

hv.extension('bokeh')
home = os.environ["HOME"]
rc("text", usetex=False)

from dask.distributed import Client
Client()

class witParticle(JITParticle):
    alphabar = Variable('alphabar', dtype=np.float32, initial=0.)
    k = Variable('k', dtype=np.float32, initial=0.)
    H = Variable('H', dtype=np.float32, initial=0.)
    vk = Variable('vk', dtype=np.float32, initial=0.)
    dvdn = Variable('dvdn',  dtype=np.float32, initial=0.)
    dkds = Variable('dkds', dtype=np.float32, initial=0.)
    dkds_dia = Variable('dkds_dia', dtype=np.float32, initial=0.)
    
    div_topo = Variable('div_topo', dtype=np.float32, initial=0.)
    div_rotary = Variable('div_rotary', dtype=np.float32, initial=0.)
    
    curv_hadv = Variable("curv_hadv", dtype = np.float32, initial = 0)
    curv_cor = Variable("curv_cor", dtype = np.float32, initial = 0)
    curv_drag = Variable("curv_drag", dtype = np.float32, initial = 0)
    curv_visc = Variable("curv_visc", dtype = np.float32, initial = 0)
    curv_prsgrd = Variable("curv_prsgrd", dtype = np.float32, initial = 0)
    curv_rate = Variable("curv_rate", dtype = np.float32, initial = 0)
    curv_total = Variable("curv_total", dtype = np.float32, initial = 0)
    
    curv_drag_sltq = Variable("curv_drag_sltq", dtype = np.float32, initial = 0)
    curv_drag_diss = Variable("curv_drag_diss", dtype = np.float32, initial = 0)
    curv_drag_sptq = Variable("curv_drag_sptq", dtype = np.float32, initial = 0)
    
    curv_adv_stretch = Variable("curv_adv_stretch", dtype = np.float32, initial = 0)
    curv_adv_divshear = Variable("curv_adv_divshear", dtype = np.float32, initial = 0)
    
    
def Sample(particle, fieldset, time):
    particle.alphabar = fieldset.alphabar[time, particle.depth, particle.lat, particle.lon]
    particle.H = fieldset.H[time, particle.depth, particle.lat, particle.lon]
    particle.k = fieldset.k[time, particle.depth, particle.lat, particle.lon]
    particle.vk = fieldset.vk[time, particle.depth, particle.lat, particle.lon]
    particle.dvdn = fieldset.dvdn[time, particle.depth, particle.lat, particle.lon]
    particle.dkds = fieldset.dkds[time, particle.depth, particle.lat, particle.lon]
    particle.dkds_dia = fieldset.dkds_dia[time, particle.depth, particle.lat, particle.lon]
    
    particle.div_topo = fieldset.div_topo[time, particle.depth, particle.lat, particle.lon]
    particle.div_rotary = fieldset.div_rotary[time, particle.depth, particle.lat, particle.lon]
    
    particle.curv_hadv = fieldset.curv_hadv[time, particle.depth, particle.lat, particle.lon]
    particle.curv_cor = fieldset.curv_cor[time, particle.depth, particle.lat, particle.lon]
    particle.curv_drag = fieldset.curv_drag[time, particle.depth, particle.lat, particle.lon]
    particle.curv_visc = fieldset.curv_visc[time, particle.depth, particle.lat, particle.lon]
    particle.curv_prsgrd = fieldset.curv_prsgrd[time, particle.depth, particle.lat, particle.lon]
    particle.curv_rate = fieldset.curv_rate[time, particle.depth, particle.lat, particle.lon]
    particle.curv_total = fieldset.curv_total[time, particle.depth, particle.lat, particle.lon]
    
    particle.curv_drag_sltq = fieldset.curv_drag_sltq[time, particle.depth, particle.lat, particle.lon]
    particle.curv_drag_diss = fieldset.curv_drag_diss[time, particle.depth, particle.lat, particle.lon]
    particle.curv_drag_sptq = fieldset.curv_drag_sptq[time, particle.depth, particle.lat, particle.lon]
    
    particle.curv_adv_stretch = fieldset.curv_adv_stretch[time, particle.depth, particle.lat, particle.lon]
    particle.curv_adv_divshear = fieldset.curv_adv_divshear[time, particle.depth, particle.lat, particle.lon]

def DeleteParticle(particle, fieldset, time):
    print("Deleting particle")
    particle.delete()

recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
            ErrorCode.ErrorThroughSurface: DeleteParticle}

state_terms = ["alphabar", "H",  "k" , "vk", "dvdn", "dkds", "dkds_dia"] 
div_terms = ["div_topo", "div_rotary"]
curv_terms = ["curv_hadv","curv_cor","curv_drag","curv_visc","curv_prsgrd","curv_rate","curv_total"]
drag_terms = ["curv_drag_sltq","curv_drag_diss","curv_drag_sptq"]
adv_terms = ["curv_adv_stretch", "curv_adv_divshear"]
variables = state_terms + curv_terms + div_terms + drag_terms + adv_terms

def interp2path(dm, theta, pathfile, r = 12.05e3):

    dm["theta_opt"] = theta
    x0, y0 = r*np.cos(theta), r*np.sin(theta)

    x,y, = dm.x_psi, dm.y_psi
    mesh = "flat"
    velocities = {'U': 'ubar_eastward',
                  'V': 'vbar_northward'}

    dimensions = {'U': {'lon': 'x_psi', 'lat': 'y_psi'},
                  'V': {'lon': 'x_psi', 'lat': 'y_psi'}}

    fieldset = FieldSet.from_xarray_dataset(dm, variables = velocities, dimensions = dimensions, mesh = mesh, time_periodic = False)

    for v in variables:
        field = Field(name = v, data = dm[v].values, lon = x, lat = y, transpose = False, mesh = mesh, allow_time_extrapolation = True)
        fieldset.add_field(field)

    pset = ParticleSet.from_list(fieldset = fieldset, pclass = witParticle, time = dm.ocean_time.values, lon = x0, lat = y0 )
    output_file = pset.ParticleFile(name = "output.nc", outputdt=delta(seconds = 300))
    kernels = AdvectionRK45 + pset.Kernel(Sample)
    pset.execute(kernels, runtime = delta(hours = 60.0), dt = delta(seconds = 300), output_file = output_file, recovery = recovery)
    output_file.export()
    dp = xr.open_dataset("output.nc")
    dp = make_dataset(dp)
    dp.to_netcdf(pathfile)
    dp.close()
    return

def make_dataset(dp):
    #dp = pd.read_pickle(filename)
    # dp["s"] = np.cumsum(np.hypot(np.gradient(dp.lon, axis = 1), np.gradient(dp.lat, axis = 1)))
    dp["s"] = np.cumsum(np.hypot(dp.lon.differentiate("obs"),dp.lat.differentiate("obs")))
    dp["V"] = dp.vk/dp.k
    dp["dhds"] = -dp.div_topo*dp.H/dp.V
    dp["alphabar"]*= (180/np.pi)
    dp["curv_drag_compute"] = (dp.curv_drag_sltq + dp.curv_drag_sptq + dp.curv_drag_diss)

    dp["alpha_path"] = np.arctan2(dp.lat.differentiate("obs"),dp.lon.differentiate("obs")) + np.pi/2
    dp["dads_path"] = dp.alpha_path.differentiate("obs")/dp.s.differentiate("obs")
    dp["d2vdns_path"] = dp.dvdn.differentiate("obs")/dp.s.differentiate("obs")
    dp["dkds_path"] = dp.k.differentiate("obs")/dp.s.differentiate("obs")
    dp["curv_adv_compute"] = (-dp.dkds_path + dp.curv_adv_divshear + dp.curv_adv_stretch)
    return dp

def domain(ds):
    x = ds.x_psi.values
    y = ds.y_psi.values
    xo = np.concatenate([x[0,:-1], x[:-1,-1], x[-1,::-1], x[-2:0:-1,0]])
    yo = np.concatenate([y[0,:-1], y[:-1,-1], y[-1,::-1], y[-2:0:-1,0]])
    coords = np.vstack((xo,yo)).T
    return hv.Curve(coords).opts(color = "k", line_width = 1)

def create_figures(df, dm):    
    ps = df.hvplot.points(x = "lon", y = "lat", s = 1, color = "r")
    pc = dm.V.hvplot.quadmesh("x_psi", "y_psi", cmap = cmocean.cm.gray_r, rasterize = True, colorbar =  False, clim = (0,.1)).opts(aspect = "equal")
    fig_xy = pc*ps*domain(dm)

    plot_dia = ( df.hvplot.line(x = "s", y = ["dkds","dkds_dia"],line_dash = ["solid","dashed"], color = ["k","k"], line_alpha = [1,1])*\
                 df.hvplot.line(x = "s", y = curv_terms, ylim = (-2e-6, 2e-6))).opts(title = "Curvature Evolution") 
    
    plot_drag = (df.hvplot.line(x = "s", y = drag_terms )*
                 df.hvplot.line(x = "s", y = "curv_drag_compute", color = "black", line_dash = "dashed", label = "curv_drag_compute" )*
                 df.hvplot.line(x = "s", y = "curv_drag", color = "black", label = "curv_drag")*
                 hv.HLine(0).opts(line_width = 1, color = "k", alpha = .25)).opts(title = "Drag Components", ylim = (-2e-6,2e-6))
  
    plot_hadv = (df.hvplot.line(x = "s", y = adv_terms + ["dkds_path"] )*
                 df.hvplot.line(x = "s", y = "curv_adv_compute", color = "black", line_dash = "dashed", label = "curv_adv_compute")*
                 df.hvplot.line(x = "s", y = "curv_hadv", color = "black", label = "curv_adv")*
                 hv.HLine(0).opts(line_width = 1, color = "k", alpha = .25)).opts(title = "Nonlinear Components", ylim = (-2e-6,3e-6))
                    
    # fig_s = (
    # plot_dia
    # +
    # df.hvplot.line(x = "s", y = ["dvdn","vk"], label ='Vorticity')
    # +
    # (df.hvplot.line(x = "s", y = drag_terms )*
    # df.hvplot.line(x="s", y = "curv_drag_compute", color = "black", line_dash = "dashed", label = "curv_drag_compute" )*
    # df.hvplot.line(x = "s", y = "curv_drag", color = "black", label = "curv_drag")*
    # hv.HLine(0).opts(line_width = 1, color = "k", alpha = .25)).opts(title = "Drag Components", ylim = (-2e-6,2e-6))
    # +
    # df.hvplot.line(x = "s", y = div_terms, label ='Divergence')
    # +
    # (df.hvplot.line(x = "s", y = adv_terms + ["dkds_path"] )*
    # df.hvplot.line(x = "s", y = "curv_adv_compute", color = "black", line_dash = "dashed", label = "curv_adv_compute")*
    # df.hvplot.line(x = "s", y = "curv_hadv", color = "black", label = "curv_adv")*
    # hv.HLine(0).opts(line_width = 1, color = "k", alpha = .25)).opts(title = "Nonlinear Components", ylim = (-2e-6,3e-6))
    # +
    # df.hvplot.line(x = "s", y = "alphabar") #, label = "depth")
    # ).opts(width = 528, shared_axes=False).cols(2)
    return fig_xy, plot_dia, plot_hadv, plot_drag

