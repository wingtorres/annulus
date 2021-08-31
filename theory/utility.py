import numpy as np
import pandas as pd
from romspy.general import coriolis

def dict2df(d, rowunits = " cm", colunits = "^{{\circ}} S", transpose = False):
    
    if not transpose:
        df = pd.Series(d).unstack()
    else:
        df = pd.Series(d).unstack().transpose()
    
    for f in df.index:
        df.rename(index = {f: f"${float(f)}$ {rowunits}"}, inplace = True )
    for f in df.columns:
        df.rename(columns = {f: f"${abs(int(f))} {colunits}$"} , inplace =  True)
    return df

def tait72(hb = 0.5, gamma = 1.0):
    return hb*( (1/16)*gamma**2 + 1/(1+8/(3*gamma**2)) )

def battjes74(Hb = 0.5, gamma=1):
    return (5/16)*Hb*gamma

def stockdon06(Hs = 1, L = 75, B = 0.1):
    return 0.35*B*(Hs*L)**0.5

# def C2d(h, z0 = 0.01, kappa = .41):
#     return kappa**2 / (np.log(h/z0) - 1 )**2

def C2d(h, cd = 0.01): #, z0 = 0.01, kappa = .41):
    return cd*h

def cd_z0(z0 = 0.01, k = 0.41, zr =1):
    return k**2/np.log(zr/z0)**2

#def jet_speed(#,Hs = 0.5, Cd = 0.01, hr = 2.0, hj = 20.0, Wr = 6e3, Wj = 500, Lr = 1e3, Ll = 1e3, g = 9.81):)

class reef_model():
    
    g = 9.807

    def __init__(self,  latitude = -30, Cd = 0.0625, 𝜆 = 0.1,
                        Wj = 1e3, Wr = 6e3, Ll = 1e3, Lr = 1e3, hj = 20, hr = 2, hl = 20,
                        island_radius = 10e3):
        
        self.f0 = coriolis(latitude)
        self.Cd = Cd
        self.h0 = hl
        self.𝜆 = 𝜆 
        self.Wj = Wj
        self.Wr = Wr
        self.Ll = Ll
        self.Lr = Lr
        self.hj = hj
        self.hr = hr
        self.hl = hl
        self.Ri = island_radius + self.Ll + self.Lr
        
        self.sw = {}
        
    def depth(self, r):
        return self.h0 + self.𝜆*r
        
    def coastal_current(self, r , f = None, Ri = None):
        if f is None:
            f = self.f0
        if Ri is None:
            Ri = self.Ri
        
        return f*(r + Ri)*np.log( (r + Ri)/Ri )

    def jet_velocity(self, etar): 
        vr = (1/self.hr)*np.sqrt(self.g*etar/self.Lr)*( 
             self.Cd/self.hr**3 + 
             (1/self.Lr)*(1/self.hj)**2*(self.Wr/self.Ll - 1) +
             self.Cd*(self.Wr/self.Lr)/self.hj**3 + 
             0*(np.pi/2)*self.Wr**2/(self.hj**2 * self.Wj**3) + 
             self.Cd*self.Wr**2/(self.hj**3 * self.Wj**2) )**(-.5)
        vj = vr*(self.hr/self.hj)*(self.Wr/self.Wj)
        self.v0 = vj
        
        return vj
    
    def sw_speed(self, s, V0 = None):
        
        if V0 == None:
            V0 = self.v0
        
        V = (V0*self.h0)/(self.h0 + self.𝜆*s)
        return V
    
    def sw_speed_grad(self, s, V0 = 0.25):
        return -V0*self.𝜆/(self.h0+s*self.𝜆)

    def sw_eta(self, s, 𝜂0 = 0.0, V0 = .25):
        return    V0**2*self.h0**2*(self.Cd - self.𝜆)/(2*self.g*self.𝜆*(self.𝜆*s + self.h0)**2) \
                - V0**2*(self.Cd - self.𝜆)/(2*self.g*self.𝜆) + 𝜂0
        
    def sw_eta_slope(self, s, V0 = 0.25):
        return V0**2*self.h0**2*(self.𝜆 - self.Cd)/(self.g*(self.h0 + self.𝜆*s)**3)
    
    def dia_sw(self,dvds,V,Cd,H,dnds): # g = 9.81):
        nl = V*dvds
        drag = Cd*V**2/H
        pgf = self.g*dnds
        return nl, drag, pgf
    
    def sw_state(self, s, V0 = 0.25):
        self.sw["s"] = s
        self.sw["h"] = self.depth(s)
        self.sw["V"] = self.sw_speed(s, V0 = V0)
        self.sw["𝜂"] = self.sw_eta(s, V0 = V0)
        self.sw["dvds"] = self.sw_speed_grad(s, V0 = V0)
        self.sw["d𝜂ds"] = self.sw_eta_slope(s, V0 = V0)
        self.sw["nlr"], self.sw["drag"], self.sw["pgf"] = self.dia_sw(self.sw["dvds"], self.sw["V"], self.Cd, self.sw["h"], self.sw["d𝜂ds"])
        
#     def make_df(self, names, variables):
#         names = ["s", "depth", "speed", "𝜂", "dnds", "dvds", "nonlinear", "drag", "pgf"]
#         data = [self.s, self.H, self.V, self.𝜂, self.d𝜂ds, self.dvds, self.sw_nlr, self.sw_drg, self.sw_pgf]
#         df = pd.DataFrame.from_dict(dict(zip(names,data)))
#         return df

