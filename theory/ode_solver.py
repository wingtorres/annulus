from numba import jit, njit
import numpy as np
import pandas as pd

@njit
def orient(theta,alpha):
    """heading angle"""
    #angle = (alpha - theta)
    angle = ( (alpha-theta) + np.pi) % (2*np.pi) - np.pi     
    return angle

@njit
def spreading(k,dadn):
    return k*dadn

@njit
def nonlinear(h,dhds,k):
    """computation of nonlinear term"""
    return dhds*2*k/h

@njit
def coriolis(h,dhds,f,u):
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
def curvature_evolution(s, Y, 𝛽 = 0.01, f = -1e-4, Cd = 2.5e-3, ri = 12e3, dadn = 0): 
    """system of ODEs"""
    𝛼, r, 𝜃, h, u, k = Y
    angle = orient(𝜃,𝛼)
    c = np.exp(1j*angle)
    drds, drdn = c.real, c.imag
    dhds = 𝛽*drds
    
    return [k, drds, (1/r)*drdn, dhds, -u*(dhds/h + dadn), 
            spreading(k,dadn) + nonlinear(h,dhds,k) + coriolis(h,dhds,f,u) + slope_torque(Cd,-𝛽*drdn,h) + dissipation(Cd,k,h)]

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
        self.S["dhds"] = self.𝛽*np.cos(self.S["angle"]) 
        self.S["dhdn"] = -self.𝛽*np.sin(self.S["angle"])
        
        self.S = pd.DataFrame.from_dict(self.S, dtype = np.float64())
        
        return self.S
    
    def diagnostics(self):
        """recover curvature equation values from ODE solution"""
        
        self.D["s"] = self.S["s"]
        
        self.D["spreading"] = spreading(self.S["k"].values, self.dadn)
        self.D["nonlinear"] = nonlinear(self.S["h"].values, self.S["dhds"].values, self.S["k"].values)
        self.D["coriolis"] = coriolis(self.S["h"].values, self.S["dhds"].values, self.f, self.S["u"].values)
        self.D["slope_torque"] = slope_torque(self.Cd, self.S["dhdn"].values, self.S["h"].values)
        self.D["dissipation"] = dissipation(self.Cd, self.S["k"].values, self.S["h"].values)
        self.D["dkds_diagnostic"] = self.D["spreading"] + self.D["coriolis"] + self.D["nonlinear"] + self.D["slope_torque"] + self.D["dissipation"]
        self.D = pd.DataFrame.from_dict(self.D, dtype = np.float64())
        
        return self.D
