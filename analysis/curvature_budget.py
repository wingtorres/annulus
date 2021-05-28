import os
import math
import numpy as np

import holoviews as hv
import hvplot.xarray
import hvplot.pandas
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import rc
import cmocean
from operator import attrgetter

from parcels import AdvectionRK4, AdvectionRK45, ErrorCode, Variable, Field, FieldSet, JITParticle, ParticleFile, ParticleSet
from datetime import timedelta as delta

import pycoawst.tools.grid as pcg
import pycoawst.tools.circulation as pcc
import pycoawst.tools.momentum as pcm

hv.extension('bokeh')
home = os.environ["HOME"]
rc("text", usetex=False)

#use setattr for witParticle and Sample

class witParticle(JITParticle):
    s = Variable('s', dtype = np.float32, initial=0.)
    lon_lag = Variable("lon_lag", dtype = np.float32, to_write = False, initial = attrgetter("lon"))
    lat_lag = Variable("lat_lag", dtype = np.float32, to_write = False, initial = attrgetter("lat"))

    alphabar = Variable('alphabar', dtype=np.float32, initial=attrgetter("alphabar"))
    k = Variable('k', dtype=np.float32, initial=attrgetter("k"))
    H = Variable('H', dtype=np.float32, initial=attrgetter("H"))
    v = Variable('v', dtype=np.float32, initial=attrgetter("V"))
    ω = Variable('ω', dtype=np.float32, initial=attrgetter("ω"))

    dadn = Variable('dadn', dtype=np.float32, initial=attrgetter("dadn"))
    dvdn = Variable('dvdn',  dtype=np.float32, initial=attrgetter("dvdn"))
    dkds = Variable('dkds', dtype=np.float32, initial=attrgetter("dkds"))
    dkds_dia = Variable('dkds_dia', dtype=np.float32, initial=attrgetter("dkds_dia"))
    
    divergence = Variable('divergence', dtype=np.float32, initial=attrgetter("divergence"))
    div_hadv = Variable('div_hadv', dtype=np.float32, initial=attrgetter("div_hadv"))
    #div_topo = Variable('div_topo', dtype=np.float32, initial=attrgetter("div_topo"))
    #div_rotary = Variable('div_rotary', dtype=np.float32, initial=attrgetter("div_rotary"))
    
    curv_hadv = Variable("curv_hadv", dtype = np.float32, initial = attrgetter("curv_hadv"))
    curv_cor = Variable("curv_cor", dtype = np.float32, initial = attrgetter("curv_cor"))
    curv_drag = Variable("curv_drag", dtype = np.float32, initial = attrgetter("curv_drag"))
    curv_visc = Variable("curv_visc", dtype = np.float32, initial = attrgetter("curv_visc"))
    curv_prsgrd = Variable("curv_prsgrd", dtype = np.float32, initial = attrgetter("curv_prsgrd"))
    curv_rate = Variable("curv_rate", dtype = np.float32, initial = attrgetter("curv_rate"))
    curv_total = Variable("curv_total", dtype = np.float32, initial = attrgetter("curv_total"))
    
    curv_drag_sltq = Variable("curv_drag_sltq", dtype = np.float32, initial = attrgetter("curv_drag_sltq"))
    curv_drag_diss = Variable("curv_drag_diss", dtype = np.float32, initial = attrgetter("curv_drag_diss"))
    curv_drag_sptq = Variable("curv_drag_sptq", dtype = np.float32, initial = attrgetter("curv_drag_sptq"))
    
    curv_adv_sheardiv = Variable("curv_adv_sheardiv", dtype = np.float32, initial = attrgetter("curv_adv_sheardiv"))
    curv_adv_veer = Variable("curv_adv_veer", dtype = np.float32, initial = attrgetter("curv_adv_veer"))
    curv_adv_shear = Variable("curv_adv_shear", dtype = np.float32, initial = attrgetter("curv_adv_shear"))
    curv_adv_curv = Variable("curv_adv_curv", dtype = np.float32, initial = attrgetter("curv_adv_curv"))
    
def Sample(particle, fieldset, time):
    dx = particle.lon - particle.lon_lag
    dy = particle.lat - particle.lat_lag
    ds = math.sqrt( math.pow(dx, 2) + math.pow(dy, 2) )
    particle.s += ds
    particle.lon_lag = particle.lon
    particle.lat_lag = particle.lat

    particle.alphabar = fieldset.alphabar[time, particle.depth, particle.lat, particle.lon]
    particle.H = fieldset.H_psi[time, particle.depth, particle.lat, particle.lon]
    particle.k = fieldset.k[time, particle.depth, particle.lat, particle.lon]
    particle.v = fieldset.v[time, particle.depth, particle.lat, particle.lon]
    particle.ω = fieldset.ω[time, particle.depth, particle.lat, particle.lon]
    particle.dvdn = fieldset.dvdn[time, particle.depth, particle.lat, particle.lon]
    particle.dadn = fieldset.dadn[time, particle.depth, particle.lat, particle.lon]
    #particle.dhds = fieldset.dhds[time, particle.depth, particle.lat, particle.lon]
    particle.dkds = fieldset.dkds[time, particle.depth, particle.lat, particle.lon]
    particle.dkds_dia = fieldset.dkds_dia[time, particle.depth, particle.lat, particle.lon]
    
    particle.divergence = fieldset.divergence[time, particle.depth, particle.lat, particle.lon]
    particle.div_hadv = fieldset.div_hadv[time, particle.depth, particle.lat, particle.lon]
    #particle.div_topo = fieldset.div_topo[time, particle.depth, particle.lat, particle.lon]
    #particle.div_rotary = fieldset.div_rotary[time, particle.depth, particle.lat, particle.lon]
    
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
    
    particle.curv_adv_shear = fieldset.curv_adv_shear[time, particle.depth, particle.lat, particle.lon]
    particle.curv_adv_sheardiv = fieldset.curv_adv_sheardiv[time, particle.depth, particle.lat, particle.lon]
    particle.curv_adv_veer = fieldset.curv_adv_veer[time, particle.depth, particle.lat, particle.lon]
    particle.curv_adv_curv = fieldset.curv_adv_curv[time, particle.depth, particle.lat, particle.lon]
    
    if particle.s >= 4e3:
        particle.state = 4
#       particle.state = 4
#       print("Halting execution...")

# def DeleteParticle(particle, fieldset, time):
#     print("Deleting particle")
#     particle.delete()

def DeleteParticle(particle, fieldset, time):
    particle.state = 4
    #print("Halting Execution: Out of bounds")

recovery = {ErrorCode.ErrorOutOfBounds: DeleteParticle,
            ErrorCode.ErrorThroughSurface: DeleteParticle}

state_terms = ["alphabar", "v", "H_psi", "k", "ω", "dvdn", "dadn", "dkds", "dkds_dia"] 
div_terms = ["divergence", "div_hadv"] #["div_topo", "div_rotary"]
curv_terms = ["curv_hadv","curv_cor","curv_drag","curv_visc","curv_prsgrd","curv_rate","curv_total"]
drag_terms = ["curv_drag_sltq","curv_drag_diss","curv_drag_sptq"]
adv_terms = ["curv_adv_shear", "curv_adv_sheardiv","curv_adv_veer", "curv_adv_curv"]
variables = state_terms + curv_terms + div_terms + drag_terms + adv_terms

def interp2path(dm, theta, pathfile, r = 12.025e3):
    dm["v"] = dm.V
    dm["theta_opt"] = theta
    x0, y0 = r*np.cos(theta), r*np.sin(theta)

    x,y, = dm.x_psi, dm.y_psi
    mesh = "flat"
    velocities = {'U': 'ubar_eastward',
                  'V': 'vbar_northward'}

    dimensions = {'U': {'lon': 'x_psi', 'lat': 'y_psi'},
                  'V': {'lon': 'x_psi', 'lat': 'y_psi'}}

    fieldset = FieldSet.from_xarray_dataset(dm, variables = velocities, dimensions = dimensions, mesh = mesh, time_periodic = False)

    for var in variables:
        field = Field(name = var, data = dm[var].values, lon = x, lat = y, transpose = False, mesh = mesh, allow_time_extrapolation = True, interp_method = "linear")
        fieldset.add_field(field)

    pset = ParticleSet.from_list(fieldset = fieldset, pclass = witParticle, time = dm.ocean_time.values, lon = x0, lat = y0 )
    output_file = pset.ParticleFile(name = pathfile, outputdt = delta(seconds = 300), convert_at_end = True)
    kernels = AdvectionRK4 + pset.Kernel(Sample)
    pset.execute(kernels, runtime = delta(hours = 120.0), dt = delta(seconds = 30), output_file = output_file, recovery = recovery, verbose_progress = True)
    display(pset)
    output_file.export()
    #output_file.close()
    return

def make_dataset(dp):
    #dp = pd.read_pickle(filename)
    # dp["s"] = np.cumsum(np.hypot(np.gradient(dp.lon, axis = 1), np.gradient(dp.lat, axis = 1)))
    #dp["s"] = np.cumsum(np.hypot(dp.lon.differentiate("obs"),dp.lat.differentiate("obs")))
    #dp["V"] = dp.vk/dp.k
    #dp["dhds"] = -dp.div_topo*dp.H/dp.V
    dp["V"] = dp.v
    dp["vk"] = dp.v*dp.k
    dp["dhds"] = dp.H.differentiate("obs")/dp.s.differentiate("obs")

    dp["alphabar"]*= (180/np.pi)
    dp["div_topo"] = -dp.v*dp.dhds/dp.H
    dp["div_rotary"] = -dp.v*dp.dadn

    dp["theta_path"] = np.arctan2(dp.lat.differentiate("obs"),dp.lon.differentiate("obs")) #+ np.pi/2
    dp["dads_path"] = dp.alpha_path.differentiate("obs")/dp.s.differentiate("obs")
    dp["sheardiv"] = dp.dvdn.differentiate("obs")/dp.s.differentiate("obs")
    dp["dkds_path"] = dp.k.differentiate("obs")/dp.s.differentiate("obs")
    
    #dp["curv_adv_compute"] = (-dp.dkds_path + dp.curv_adv_sheardiv + dp.curv_adv_stretch)
    #dp["curv_drag_compute"] = (dp.curv_drag_sltq + dp.curv_drag_sptq + dp.curv_drag_diss)
    return dp

def create_figures(df, dm):    
    #ps = df.hvplot.points(x = "lon", y = "lat", s = 1, color = "r")
    #pc = dm.V.hvplot.quadmesh("x_psi", "y_psi", cmap = cmocean.cm.gray_r, rasterize = True, colorbar =  False, clim = (0,.1)).opts(aspect = "equal")
    #fig_xy = pc*ps*domain(dm)

    plot_vrt =  df.hvplot.line(x = "s", y = ["vk","dvdn"],  color = ["orange","purple"]).opts(title = "Vorticity")#*hv.Hline(0).opts(line_width = 1, color = "k", alpha = .25) 

    # plot_div =  
    

    

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
    return plot_dia, plot_hadv, plot_drag, plot_div, plot_vrt

