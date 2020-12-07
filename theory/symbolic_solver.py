from sympy import Function, Eq, Derivative, symbols, solve, lambdify, plot, pi, powsimp, latex, init_printing
from sympy.solvers.ode import dsolve, classify_ode
from sympy.abc import r

r = symbols('r', cls=Function)
s = symbols('s')
R, U, f, C, B = symbols('R U f C_D B', real = True, constant = True)

def analytic_solution():

	expr = Eq(r(s).diff(s,3), r(s).diff(s,1)*(1/R**2 + f/(U*R)) - 2*C/(B*R**2))
	sol_par = dsolve(expr, r(s))

	ics = [sol_par.rhs.subs(s,0) - R, 
		   sol_par.rhs.diff(s,1).subs(s,0) + pi/2, 
		   sol_par.rhs.diff(s,2).subs(s,0) - 1/R + pi/2/R**2]
	C1,C2,C3 = symbols('C1 C2 C3')
	coeffs = solve(ics,[C1,C2,C3])

	sol = sol_par.subs(coeffs)
	sol_exa = sol.simplify()
	sol_exa_disp = powsimp(sol, deep = True, force = False, combine = 'exp')
	return  latex(expr), latex(sol_par), latex(sol_exa_disp), sol_exa

def analytic_function(sol,radius,coriolis,speed,slope,drag = 0):
	return	lambdify(s, sol.rhs.subs(R,radius).subs(f,coriolis).subs(U,speed).subs(B,slope).subs(C,drag),'numpy')
