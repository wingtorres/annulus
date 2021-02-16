from numba import jit, njit
import numpy as np
import pandas as pd
import xarray as xr
#import pycoawst.tools.circulation as pcc

@njit
def orient(theta,alpha):
    """heading angle"""
    #angle = (alpha - theta)
    angle = ( (alpha-theta) + np.pi) % (2*np.pi) - np.pi     
    return angle

#https://stackoverflow.com/questions/25829364/applying-a-decorator-to-an-imported-function
# nonlinear_topo = njit(pcc.nonlinear_topo)
# nonlinear_spr = njit(pcc.nonlinear_spr)
# coriolis_topo = njit(pcc.coriolis_topo)
# coriolis_spr = njit(pcc.coriolis_spr)
# slope_torque = njit(pcc.slope_torque)
# dissipation = njit(pcc.dissipation) 

@njit
def spreading(k,dadn):
    return k*dadn

@njit
def nonlinear_stretching(h,dhds,k):
    """computation of nonlinear term"""
    return dhds*2*k/h

@njit
def nonlinear_spreading(dadn,k):
    """computation of nonlinear term"""
    return k*dadn

@njit
def coriolis_stretching(h,dhds,f,u):
    """computation of coriolis term"""
    return dhds*f/(u*h)

@njit
def slope_torque(Cd,dhdn,h):
    """computation of slope torque term"""
    return -Cd*dhdn/(h*h)

@njit
def dissipation(Cd,k,h):
    """computation of dissipation term"""
    return -Cd*k/h

@njit
def stokes_divergence(h, dhds, Hs = 2.0, k = 2*np.pi/100, omega = 2*np.pi/10):
    """
    Depth-averaged stokes drift velocity
    
    "Stokes transport in layers in the water column based
    on long-term wind statistics", Myrhaug et al. 2018 Oceanologia

    """
    #us = 9.81*(Hs/2)**2*k/(2*omega*h) #*(np.sinh(2*k*d) / np.sinh(2*k*d))

    return Hs**2*k/(2*omega*h**2)*dhds

@njit
def curvature_evolution(s, Y, 𝛽 = 0.01, f = -1e-4, Cd = 2.5e-3, ri = 12e3, dadn = 0): 
    """system of ODEs"""
    𝛼, r, 𝜃, h, u, k = Y
    angle = orient(𝜃,𝛼)
    c = np.exp(1j*angle)
    drds, drdn = c.real, c.imag
    dhds = 𝛽*drds
    dhdn = -𝛽*drdn
    dkds = (
            nonlinear_stretching(h,dhds,k) + 
            nonlinear_spreading(dadn,k) + 
            coriolis_stretching(h,dhds,f,u) + 
            slope_torque(Cd,dhdn,h) + 
            dissipation(Cd,k,h) 
            # + (2*k*u + f)*stokes_divergence(h, dhds)/u**2
            )/5

    return [k, drds, (1/r)*drdn, dhds, -u*dhds/h, dkds]

@njit
def reefcrest(s,Y,𝛽,f,Cd,ri,dadn): 
    """stop integration if jet trajectory runs into reefcrest"""
    return (Y[1] - ri  + 1)
reefcrest.terminal = True

@njit
def leavedomain(s,Y,𝛽,f,Cd,ri,dadn): 
    """stop integration if jet trajectory runs into reefcrest"""
    return (-(Y[1] - 15e3  + 1))
leavedomain.terminal = True

class ODEsol():
    """custom class to store parameters and solutions for ODE solutions w/ methods 
    to store the state variables and calculate diagnostics or other derived quantities
    """
    
    def __init__(self, 𝛽, f, Cd, ri = 12e3, dadn = 0):
        """input parameters"""
        
        self.𝛽 = 𝛽
        self.f = f
        self.Cd = Cd
        self.dadn = dadn
        self.ri = 12e3
        self.h0 = 20
        self.start = 0
        self.finis = 25e3
        self.S = {}
        self.K = {}
        self.D = {}
        
    def state(self,sol,s):
        """load state variables"""
        
        sol = sol.sol(s)
        self.s = s
        self.S["s"] = s#/1e3
        self.S["𝛼"], self.S["r"], self.S["𝜃"]  = sol[0,:], sol[1,:]/1e3, sol[2,:]
        self.S["x"], self.S["y"] = self.S["r"]*np.cos(self.S["𝜃"]), self.S["r"]*np.sin(self.S["𝜃"])
        self.S["h"], self.S["u"], self.S["k"] = sol[3,:], sol[4,:], sol[5,:]
        self.S["angle"] = orient(self.S["𝜃"] , self.S["𝛼"])
        self.S["ω"] = self.S["u"]*self.S["k"]
        self.S["dhds"] = self.𝛽*np.cos(self.S["angle"]) 
        self.S["dhdn"] = -self.𝛽*np.sin(self.S["angle"])
        
        self.S = pd.DataFrame.from_dict(self.S, dtype = np.float64())
        
        return self.S

    def diagnostics(self):
        """recover curvature equation values from ODE solution"""
        
        s = self.S["s"]
        #https://github.com/pydata/xarray/issues/2931
        dp = xr.Dataset(data_vars = {"x": ('s', self.S["x"]), "y": ('s', self.S["y"])}, coords = {'s': ('s',s)}).set_coords(['x','y'])
            
        dp["nonlinear_stretching"] = ('s', nonlinear_stretching(self.S["h"].values, self.S["dhds"].values, self.S["k"].values) )
        dp["nonlinear_spreading"] = ('s', nonlinear_spreading(self.dadn, self.S["k"].values) )
        dp["coriolis_stretching"] = ('s', coriolis_stretching(self.S["h"].values, self.S["dhds"].values, self.f, self.S["u"].values) )
        dp["slope_torque"] = ('s', slope_torque(self.Cd, self.S["dhdn"].values, self.S["h"].values) )
        dp["dissipation"] = ('s', dissipation(self.Cd, self.S["k"].values, self.S["h"].values) )

        # self.D["s"] = self.S["s"]
        # self.D["nonlinear_topo"] = nonlinear_topo(self.S["h"].values, self.S["dhds"].values, self.S["k"].values)
        # self.D["nonlinear_spr"] = nonlinear_spr(self.dadn, self.S["k"].values)we
        # self.D["coriolis_topo"] = coriolis_topo(self.S["h"].values, self.S["dhds"].values, self.f, self.S["u"].values)
        # self.D["coriolis_spr"] = coriolis_spr(self.dadn, self.f, self.S["u"].values)
        # self.D["slope_torque"] = slope_torque(self.Cd, self.S["dhdn"].values, self.S["h"].values)
        # self.D["dissipation"] = dissipation(self.Cd, self.S["k"].values, self.S["h"].values)

        dp["dkds_diagnostic"] = dp.nonlinear_stretching + dp.nonlinear_spreading + dp.coriolis_stretching + dp.slope_torque + dp.dissipation
        #dp["dkds_diagnostic"] = dp.nonlinear_topo + dp.nonlinear_spr + dp.coriolis_topo + dp.coriolis_spr + dp.slope_torque + dp.dissipation
        #self.D["dkds_diagnostic"] = self.D["spreading"] + self.D["coriolis"] + self.D["nonlinear"] + self.D["slope_torque"] + self.D["dissipation"]
        self.D = dp.to_dataframe() #pd.DataFrame.from_dict(self.D, dtype = np.float64())
        
        return self.D
