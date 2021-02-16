import panel as pn
import numpy as np

navbar = """
[Home](/home) | 
[2D Numerical Solution](/diagnostics_2D) | 
[ODE Numerical Solution](/ode_solution-app) | 
[ODE Stability](/ode_stability-app)
"""

navpane = pn.pane.Markdown(navbar)

def m2km(value):
    return f"{abs(value)/1e3:d} km"

def border(ds):
    x = ds.x_psi.values
    y = ds.y_psi.values
    xo = np.concatenate([x[0,:-1], x[:-1,-1], x[-1,::-1], x[-2:0:-1,0]])
    yo = np.concatenate([y[0,:-1], y[:-1,-1], y[-1,::-1], y[-2:0:-1,0]])
    coords = np.vstack((xo,yo)).T
    #hv.Curve(coords).opts(color = "k", line_width = 1)
    return coords

