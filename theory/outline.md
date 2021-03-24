# Curvature Dynamics of a Coastal Barotropic Jet on a Slope

## Abstract

## 1.0 Introduction

- **Outflow onto a slope is a common geophysical phenomenon encompassing a variety of flows including rivers/streams, engineering discharges, tidal jets, rip currents, and wave-driven flows on reefs** (Horner-Devine et al. 2015, Jirka & Donekeer 1991, Wolanski et al. 1988, Symonds et al. 1995). These circulation patterns are characterized by the presence of a steep gradient in water properties such as momentum and density relative to the ambient fluid, and generally play a key role in the transport of scalars in the coastal ocean. As the scale of these features

Another feature of these flows is that the governing dynamics tend to change along their course, typically understood as a transition nearfield to farfield dynamics as the length/time scales develop. The nearfield tends to be inertia-dominated, while the farfield tends to be rotationally-dominated; often the circulation pattern evolves into geostrophic balance as an alongshore coastal current. The mechanics of this transition are often complicated by the overlapping effects of buoyancy, bathymetry, alongshore currents, tides, winds, and mixing.

- **This study focuses on developing theory regarding the trajectory of consecutive small-intermediate scale barotropic jets in shallow water, recognizing the considerable emphasis that has been placed on large-scale single outflows, especially those subject to stratification, such as river plumes**. Weakly stratified small-scale features such as tidal jets and reef pass jets are physically and ecologically important, but have not often been approached from a geophysical perspective as sometimes the physical scaling arguments do not justify accounting for the Coriolis acceleration (Wolanski et al. 1988 CSR, Herdman et al. 2016 JGR). Here, however, we identify realistic conditions for which rotation may play a non-trivial role in the behavior of outflow jets at scales far smaller than the Rossby radius of deformation. These features tend to be geophysically "thin" - the jet width is much smaller than the radius of curvature and the inertial radius, and so we justify taking inspiration from select larger-scale GFD studies in order to explain small-scale dynamics. Interestingly, these conditions lie within the scale gap between geophysical and engineering literature as identified by Basdurak 2020 JGR - the former generally considering wide aspect $\mathcal{O}(0.001-0.01)$ rotationally-dominated outflows and the latter generally considering narrow aspect $\mathcal{O}(0.1-1)$ inertia-dominated outflows. Based on scaling arguments and in the abscence of stratification we characterize these features using the Rossby number and equivalent Reynold's number which depend on the bottom friction, rotation, and nonlinear advection forces. 

- **Foundational work on the path of free jets (not necessarily outflows) by Warren 1963, interested in Gulf Stream meanders, demonstrated how the application of the barotropic vorticity equation could be used to derive a path equation for the trajectory of a free jet in the presence of topography**. This was later extended to arbitrary topography by Cushman-Roisin 1997 who also characterized new behaviors such as jet trapping, retroflection, and deflection. Key to their efforts was solving for curvature in the vorticity equation, which arises when cast in streamwise normal coordinates. Woods & Beardsley 1988 CSR numerically simulated a 2D barotropic outflow onto a linear slope (the subject at hand here) and identified how nonlinear topographic stretching results in nearfield acquisition of cyclonic vorticity, first deflecting the jet opposite to its inertial tendency. Eventually bottom friction damps this nonlinearity and the jet develops into a downcoast current, but the relevant contribution of the work was to the importance of nonlinear vorticity dynamics in the nearfield. In this manuscript we synthesize these approaches to the problem, using curvature via the vorticity equation to understand the dynamics of 2D barotropic flow on a slope.
- Research Questions
    - How do slope/bottom friction/Coriolis effect
interact to influence the trajectory of coastal outflows?
    - Can vorticity/curvature dynamics explain the evolution of the jet developing over topography?
    - Under what conditions, if it all, can the jet behavior be captured by a 1D reduced physics model?
    
- **Here, we investigate the behavior of a barotropic coastal jet outflow by deriving a novel equation for its curvature then using it to examine an array of numerical simulations, comparing the results against 1D analytic predictions**. 
    - Theory (Section 2): Summarize existing work on vorticity equation and natural coordinates. Derive general curvature evolution equation and 1D simplification
    - Experimental design (Section 3): Describe numerical model experiments covering model configuration and parameter space
    - Results (Section 4): 
        - Detailed description of illustrative single run 
    - Summary (Section 5) findings across parameter space to identify regimes/trends. Validity of 1D simplification assessed


    


## Theory
 
### Governing Equations

The depth-integrated 2D shallow water system
in a streamwise normal coordinate system $(s,n)$

![Figure 1](https://journals.ametsoc.org/view/journals/phoc/47/5/full-jpo-d-16-0239.1-f1.jpg)
I like this diagram from Wenegrat + Thomas 2017 JPO Fig. 1. I propose we ask to reuse it with modification to accommodate our choice of variables, and perhaps add $\vec{u},\vec{v}$ vectors too?



#### Continuity
$$
\frac{\partial \eta }{\partial t} +  \frac{\partial(V H)}{\partial s} + VH\frac{\partial \alpha}{\partial n} = 0
\tag{1}
$$


#### Momentum
$$
\frac{\partial V}{\partial t} + 
     V \frac{\partial V}{\partial s} + 
    C_d\frac{V^2}{H} +
    g\frac{\partial \eta}{\partial s} = 0
\tag{2}
$$

   
$$ 
V\frac{\partial \alpha}{\partial t} + V^2 K + fV +  g\frac{\partial \eta}{\partial n} = 0 
\tag{3}
$$

where $H(s,n,t) = h(s,n) + \eta(s,n,t)$. 

#### Vorticity
Start w/ Signell + Geyer 1991 JGR

$$
\frac{\partial \omega }{\partial t} + 
\vec{u} \cdot \nabla (\omega + f) + 
(\omega + f)( \nabla \cdot \vec{u} ) + 
\frac{C_d V}{H^2}[\vec{u} \times \nabla H] -
\frac{C_d}{H}[\vec{u} \times \nabla V] +
\frac{C_d V \omega}{H} = 0 
$$

Substituting the definition of vorticity and divergence in (<i>s,n</i>) coordinates, assuming steady flow, and dividing by $V^2$ produces the curvature equation...

$$
\frac{\partial k}{\partial s} =
\frac{\partial}{\partial s}\frac{\partial V}{\partial n} +
k\frac{\partial V}{\partial s} -
\frac{1}{V}\frac{\partial f}{\partial s} +
(-\frac{\partial V}{\partial n} + VK  + f)\bigg[ \frac{\partial V}{\partial s} + V\frac{\partial \alpha}{\partial n} \bigg] -
\frac{C_d}{H^2}\frac{\partial H}{\partial n} +
\frac{C_d}{VH}\frac{\partial V}{\partial n} -
\frac{C_d}{VH}\big[-\frac{\partial V}{\partial n} + VK \big] 
\tag{4}
$$

Rearranging continuity (1) gives us an expression for $\frac{\partial V}{\partial s}$

$$ \frac{1}{V}\frac{\partial V}{\partial s} = -\frac{1}{H}\frac{\partial H}{\partial s} - \frac{\partial \alpha}{\partial n} $$

which upon substitution gives

$$
\underbrace{\frac{\partial K}{\partial s}}_{\text{curvature gradient}} =
\underbrace{ \frac{1}{V}\frac{\partial}{\partial s}\frac{\partial V}{\partial n} }_{\text{shear divergence}} + 
\underbrace{ K\frac{\partial \alpha}{\partial n} }_{\text{confluence veering}} -
\underbrace{ \frac{1}{V}\frac{\partial f}{\partial s} }_{\text{Beta effect}} +
\underbrace{\bigg[ \frac{1}{VH} \frac{\partial H}{\partial s}  \bigg]}_{\text{stretching}}\bigg[
\underbrace{-\frac{\partial V}{\partial n}}_{\text{shear}} + 
\underbrace{2VK}_{\text{curvature}} +
\underbrace{f}_{\text{Coriolis}}\bigg]  \\
-\underbrace{ \frac{C_d}{H^2}\frac{\partial H}{\partial n} }_{\text{slope torque}} +
\underbrace{ \frac{2C_d}{VH}\frac{\partial V}{\partial n} }_{\text{speed torque}} -
\underbrace{ \frac{C_d}{H}K }_{\text{dissipation}} 
$$

### Interpretation of terms
- Shear divergence: The tendency for streamlines to curve along the transition between different velocity profiles. For a momentum jet, this processes outwardly 
deflects streamlines on its flanks.
![](https://i.imgur.com/QPuqIbj.png)

- Confluence veering: The interaction of diffluent/confluent streamlines with existing flow curvature further enhances/attenuates curvature in the same direction.
![](https://i.imgur.com/mVmPRCF.png)


- Beta effect:
Curvature acquired due to the streamwise variation of the Coriolis parameter. A streamline differentially affected by rotation along its path will naturally acquire curvature.
![](https://i.imgur.com/Unhp0Tx.png)

- Nonlinear stretching 

    - Curvature: The tendency for cross-isobath flow with curvature to deflect
![](https://i.imgur.com/nr674cF.png)

    - Shear: Cross-isobath flow with shear will curve to conserve angular momentum via vortex stretching. In a momentum jet this manifests as outward lateral deflecting during it shoaling, and inward as it deepns.
![](https://i.imgur.com/re2YmZI.png)


    - Coriolis: Flow divergence due to topography implies a differential effect of the Coriolis force along streamlines inducing curvature. 
![](https://i.imgur.com/4lvWvN4.png)

- Slope Torque: In the presence of bottom drag, flow parallel to isobaths will curve shoreward due to depth-dependent bottom friction acting more strongly on the shallow side of the flow feature
![](https://i.imgur.com/ooOz00A.png)


- Speed torque:
If there is a normal gradient in speed, the flow will experience a torque due to the differential forces on either side of the streamline inducing curvature
![](https://i.imgur.com/HnnPpOJ.png)

- Dissipation:
Bottom friction will also act to reduce curvature present in the flow
![](https://i.imgur.com/8aobi3V.png)


## 1D System
We now elect to integrate equation (4) along the streamline defined by the local speed maxima along the jet's trajectory -an inflection point in the velocity profile where $\frac{\partial V}{\partial n} = 0$-consequently eliminating the terms involving velocity shear.Also considering the problem on an $\textit{f}$-plane ($\frac{\partial f}{\partial s} = 0$) gives

$$
\frac{\partial k}{\partial s} =
k\frac{\partial \alpha}{\partial n} +
(2K + f)\bigg[ \frac{1}{H} \frac{\partial H}{\partial s}  \bigg] -
\frac{C_d}{H^2}\frac{\partial H}{\partial n} +
\frac{C_d}{H}k 
$$

Finally, assuming jet streamlines are parallel along the jet trajectory $\frac{\partial \alpha}{\partial n} = 0$ leaves us with equation (5) below. Doing so reduces the dimensionality of the system of equations; state variables become a function of along trajectory distance only. Even though this assertion precludes jet spreading, we expect the emergent volume balance dynamics in the jet will be dominated by the $\frac{\partial (VH)}{\partial s}$ terms, where the outflow is fast and streamwise depth gradient steep.One choice for a non-dimensional parameter quantifying the primacy of the rotary and streamwise terms in the continuity equation is $\Upsilon = \frac{W\lambda}{2 H_0 \alpha_0}$, where $W/2$ is half of the jet width, $\lambda$ the slope, $H_0$ the initial depth and $\alpha_0$ the spreading angle. For characteristic small-scale coastal outflows over steep slopes, e.g. wave-driven flow on some Pacific atolls ($\lambda > 0.1$), $W/2 = 500 m, H_0 = 20m, \alpha_0 \sim \pi/24$ so $\Upsilon \approx 20$, whereas for weak slopes typical of beaches on continental shelves ($\lambda < 0.01$) $\Upsilon$ may be closer to unity or less.

$$
\frac{\partial k}{\partial s} =
(2K + f)\bigg[ \frac{1}{H} \frac{\partial H}{\partial s}  \bigg] -
\frac{C_d}{H^2}\frac{\partial H}{\partial n} +
\frac{C_d}{H}k 
\tag{5}
$$

For a radially symmetric skirted island, the bathymetry follows the function $h = h_0 + \lambda(r - R_i)$ where $\lambda$ is the bottom slope, $r$ is the radial coordinate relative to the island center, and $R_i$ is the island radius. Using the continuity and the curvature equations, we can now develop a system of ODE's describing flow kinematics along the streamwise coordinate of the jet.


\begin{aligned}
\frac{d\alpha}{ds} &= k \\ 
\frac{dr}{ds} &= \cos(\alpha - \theta ) \\
\frac{d\theta}{ds} &= \frac{1}{r}\sin(\alpha - \theta) \\
\frac{dh}{ds} &= \lambda\frac{dr}{ds} \\
\frac{du}{ds} &= -\frac{u}{h}\frac{dh}{ds} \\
% \frac{dk}{ds} &= \frac{f}{uh}\frac{dh}{ds} \\ 
\frac{dk}{ds} &= \frac{1}{h}\frac{dh}{ds}\big[2k + \frac{f}{u}\big] - C_D\frac{1}{h^2}\sin{(\alpha-\theta)} - C_D \frac{k}{h} \\ 
\end{aligned}

Here, α is the flow angle in cartesian coordinates while θ is the azimuthal island coordinate (Cushman-Roisin 1997 DAO). This system of equations will be integrated using the Radau IIa implicit Runge-Kutta method provided in SciPy's <i>solve_ivp</i> function.

## Experimental Design

A series of numerical experiments were carried out to explore the dynamics of the coastal outflow using the curvature equation derived above under different conditions: all 27 permutations of the following latitude $\phi$, bottom slope $\lambda$, and bottom drag coefficient $C_D$ values. The case in bold will be explored in detail.
- $\phi = [-1.0^{\circ}, -15.0^{\circ}, \mathbf{-30.0^{\circ}}]$
- $\lambda = [0.01, 0.05, \mathbf{0.1}]$
- $C_D = [0.0625, \mathbf{0.125}, 0.25]$

For each case, the terms in equation (4) will be computed along the center streamline of the jet and the leading order balances along its trajectory identified




### The Idealized Domain

The bathymetry (h) is analytically prescribed to be $\pi/6$ sector of a conical frustrum (skirted island with linearly sloping bathymetry), and so a polar grid is used.$$ h(r) = \lambda(r - r_i) + h_0 $$ where the minimum depth $h_0 = 20$m, the island radius $r_i = 12\text{km}$, and the slope $\lambda= 0.1$, typical of Pacific atolls with steep slopes. An outflow boundary condition at the inner boundary $r = r_i$ is specificed using a von Mises distribution with a negative offset 
so that the net transport across the boundary is 0: $\int^{\pi/12}_{\pi/12}v(\theta) \, d\theta = 0$. This is meant to capture how wave-driven inflow at the reef crest should match the jet outflow.

$$v(r = r_i, \theta) = V_0 \big[ \frac{e^{\kappa cos(\theta)}}{2\pi I_0(\kappa)} - b \big]$$

where $v$ is the radial velocity, $I_0$ is the zeroth order modified Bessel function, $\kappa$ is a shape parameter determining the steepness of the profile, $b$ is an offset parameter, and $\theta$ is the azimuthal coordinate. Scipy's <i>minimize</i> routine was used to find $b$ (Nelder-Mead method), while $k$ was tuned so the jet width $W_J$ = 1 km. 
![](https://i.imgur.com/KCIPz5h.png)

- Bathymetry figure?


### Numerical Model Configuration

The depth-integrated equations of motion were solved numerically on this domain with the Regional Ocean Modeling System ($\href{https://www.myroms.org/}{\text{ROMS}}$). An analytic depth-dependent quadratic bottom friction scheme was used $\tau_b = C_D\frac{V^2}{h}$ to match the governing equations used in the derivation, where $V$ is the flow speed. The harmonic lateral viscosity was set to $0.1 m^2s^{-1}$. The Coriolis parameter is constant throughout the domain (f-plane), simulations were integrated for 10 days at a timestep of $0.25s$ while the grid is spaced $15m$ in the radial direction and $15m-20m$ in the azimuthal direction. Periodic boundary conditions were imposed on the azimuthal boundaries of the domain and Flather-type boundary conditions (Shchepetkin BC's) at the radial boundaries. The inflow speed scale is $V_0 = 0.125 ms^{-1}$

### Center streamline identification

Because equation (4) is valid only along a streamline while the key simplification to arrive at (5) is $\frac{1}{V}\frac{\partial V}{\partial n} = 0, \frac{\partial \alpha}{\partial n} = 0$, the jet streamline minimizing the objective function
$$J = \frac{\int_C \big( |\frac{1}{V}\frac{\partial V}{\partial n}| + |\frac{\partial \alpha}{\partial n}| \big) \, ds}{\int_C \, ds} $$ 

is of interest. The optimal streamline was found using $\href{https://oceanparcels.org/}{\text{Ocean Parcels}}$ particle tracking package in conjunction with Scipy's <i>minimize_scalar</i> function, where the initial azimuthal location of the particle along $r_i$ was optimized. Lagrangian particles were integrated over a static velocity field i.e. the time-average over the final three hours of simulation, thus defining mean streamlines. 

## Results 

### The barotropic outflow
section focused on describing problem and displaying fields

- vorticity components figure (xy, s)
![](https://i.imgur.com/fczLlPI.png)
validity of $\frac{\partial V}{\partial n} = 0$

- divergence components figure (xy, s)
![](https://i.imgur.com/xTxr5YQ.png)
validity of $V\frac{\partial \alpha}{\partial n} = 0$

- momentum components figure (xy, s)?
- curvature components figure (xy?, s)

![](https://i.imgur.com/OiA8Yps.png)

### The effect of rotation, slope, and bottom friction on outflow dynamics
An analysis focused section explaining how the role of terms across the parameter space influences the physics of the jet. Discuss utility/effectiveness of 1D model

- Trajectories 
![](https://i.imgur.com/dg2gUKu.png)
Red is actual jet streamline and blue is 1D solution

- curvature components $s$: ($\phi, \lambda, C_D$)
![](https://i.imgur.com/XiN3XBB.png)

## Summary

## References

## Appendix


### Streamwise Normal Coordinate System

\begin{aligned}
\hat{s} = \frac{\vec{u}}{V} \\
\hat{n} = k \times \hat{k} \\
\alpha = \tan^{-1}\frac{\vec{v}}{\vec{u}} \\
\frac{\partial \alpha}{\partial s} = K \\
\nabla \cdot \vec{u} = \frac{\partial V}{\partial s} + V\frac{\partial \alpha}{\partial n} \\
\nabla \times \vec{u} = -\frac{\partial V}{\partial n} + VK \\
\end{aligned}