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
def nonlinear_topo(h,dhds,k):
    """computation of nonlinear term"""
    return dhds*2*k/h

@njit
def coriolis_topo(h,dhds,f,u):
    """computation of coriolis term"""
    return dhds*f/(u*h)

@njit
def nonlinear_spr(dadn,k):
    """computation of nonlinear term"""
    return 2*k*dadn

@njit
def coriolis_spr(dadn,f,u):
    """computation of coriolis term"""
    return dadn*f/u


@njit
def slope_torque(Cd,dhdn,h):
    """computation of slope torque term"""
    return -Cd*dhdn/(h*h)

@njit
def dissipation(Cd,k,h):
    """computation of dissipation term"""
    return -Cd*k/h

@njit
def curvature_evolution(s, Y, 𝛽 = 0.01, f = -1e-4, Cd = 2.5e-3, ri = 12e3, dadn = 0): 
    """system of ODEs"""
    𝛼, r, 𝜃, h, u, k = Y
    angle = orient(𝜃,𝛼)
    c = np.exp(1j*angle)
    drds, drdn = c.real, c.imag
    dhds = 𝛽*drds
    return [k, drds, (1/r)*drdn, dhds, -u*(dhds/h + dadn), 
            nonlinear_topo(h,dhds,k) + nonlinear_spr(dadn,k) +
            coriolis_topo(h,dhds,f,u) + coriolis_spr(dadn,f,u) + 
            slope_torque(Cd,-𝛽*drdn,h) + dissipation(Cd,k,h)]

#    return [k, drds, (1/r)*drdn, dhds, -u*(dhds/h + dadn), 
#            spreading(k,dadn) + nonlinear(h,dhds,k) + coriolis(h,dhds,f,u) + slope_torque(Cd,-𝛽*drdn,h) + dissipation(Cd,k,h)]

@njit
def reefcrest(s,Y,𝛽,f,Cd,ri,dadn): 
    """stop integration if jet trajectory runs into reefcrest"""
    return (Y[1] - ri  + 1)
reefcrest.terminal = True

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
        self.S["s"] = s/1e3
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
            
        dp["nonlinear_topo"] = ('s', nonlinear_topo(self.S["h"].values, self.S["dhds"].values, self.S["k"].values) )
        dp["nonlinear_spr"] = ('s', nonlinear_spr(self.dadn, self.S["k"].values) )
        dp["coriolis_topo"] = ('s', coriolis_topo(self.S["h"].values, self.S["dhds"].values, self.f, self.S["u"].values) )
        dp["coriolis_spr"] = ('s', coriolis_spr(self.dadn, self.f, self.S["u"].values) )
        dp["slope_torque"] = ('s', slope_torque(self.Cd, self.S["dhdn"].values, self.S["h"].values) )
        dp["dissipation"] = ('s', dissipation(self.Cd, self.S["k"].values, self.S["h"].values) )

        # self.D["s"] = self.S["s"]
        # self.D["nonlinear_topo"] = nonlinear_topo(self.S["h"].values, self.S["dhds"].values, self.S["k"].values)
        # self.D["nonlinear_spr"] = nonlinear_spr(self.dadn, self.S["k"].values)
        # self.D["coriolis_topo"] = coriolis_topo(self.S["h"].values, self.S["dhds"].values, self.f, self.S["u"].values)
        # self.D["coriolis_spr"] = coriolis_spr(self.dadn, self.f, self.S["u"].values)
        # self.D["slope_torque"] = slope_torque(self.Cd, self.S["dhdn"].values, self.S["h"].values)
        # self.D["dissipation"] = dissipation(self.Cd, self.S["k"].values, self.S["h"].values)
        dp["dkds_diagnostic"] = dp.nonlinear_topo + dp.nonlinear_spr + dp.coriolis_topo + dp.coriolis_spr + dp.slope_torque + dp.dissipation
        #self.D["dkds_diagnostic"] = self.D["spreading"] + self.D["coriolis"] + self.D["nonlinear"] + self.D["slope_torque"] + self.D["dissipation"]
        self.D = dp.to_dataframe() #pd.DataFrame.from_dict(self.D, dtype = np.float64())
        
        return self.D
