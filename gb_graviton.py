import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
import sys
from math import factorial, gamma, pi

# The final output should not have the factor 1/gtot. Instead it should have a factor 1/2. This is because blackhawk multiplies the cross sections by g = 2, the 4D graviton degrees of freedom, instead of gtot. Otherwise, all degeneracies are correctly handled.

n = int(sys.argv[1])
s = 2
rh = 1
gtot = (n+1)*(n+4)/2 # total number of degrees of freedom, equal to sum(Ns)
Ns = [1, n+1, n*(n+3)/2] # number of each type of perturbation

def f(r,n): 
	return 1 - (rh/r)**(n+1)

def rr(y,n):
	if y < -10**(-10):
		out = rh*(1-np.exp(y))**(-1/(n+1))
	else:
		out = rh*(-y)**(-1/(n+1))
	return out
r = np.vectorize(rr)

ks = [0,3,-1]
# t is type of perturbation, taking values 0, 1, or 2
def deriv(p,y,om,l,n,t): # does not yet take account of scalars
	ReR, ImR, ReS, ImS = p
	h = f(r(y,n),n)
	s = r(y,n)
	k = ks[t]
	c = (n+1)*rh**(n+1)
	V = (h/s**2)*(l*(l+n+1) + n*(n+2)/4 - k*(n+2)**2*(1-h)/4)
	if t == 0:
		m = l*(l+n+1)-n-2
		q = (n+2)**4*(n+3)**2
		z = 16*m**3+4*m**2*(n+2)*(n+4)
		pp = (n+2)*(n+3)*(4*m*(2*n**2+5*n+6)+n*(n-2)*(n+2)*(n+3))
		w = -12*m*(n+2)*(m*(n-2)+n*(n+2)*(n+3))
		Vtop = q*(1-h)**3 + pp*(1-h)**2 + w*(1-h) + z
		Vbot = 4*(2*m + (n+2)*(n+3)*(1-h))**2
		V = h*Vtop/(s**2*Vbot)
	x = s**(2*n+4)*(V - om**2)/c**2
	der = (n+2)*h*s**(n+1)/c
	return [ReS, ImS, x*ReR + der*ReS, x*ImR + der*ImS]

def degen(l,n,t): # the number of modes for each value of l (analogous to factor of 2*l+1)
	if t == 0:
		return (2*l+n+1)*factorial(l+n)/(factorial(l)*factorial(n+1))
	if t == 1:
		return l*(l+n+1)*(2*l+n+1)*factorial(l+n-1)/(factorial(l+1)*factorial(n))
	if t == 2:
		top = n*(n+3)*(l+n+2)*(l-1)*(2*l+n+1)*factorial(l+n-1)
		return top/(2*factorial(l+1)*factorial(n+1))

###########################################################
# initialising the solver
###########################################################
steps = 10**6
A = 1
ypo = []
for m in range(7):
	points = np.linspace(np.log10(5*(m+1)), -3*(m+1), steps+1)
	ypo.append(-np.power(10, points))

def grey(om,l,n,t):
	###########################################################
	# solving the ODE
	###########################################################
	omeg = om*rh/(n+1)
	p0 = [A,0,0,-A*omeg]
	ypoints = ypo[n]
	sol = odeint(deriv, p0, ypoints, args=(om,l,n,t))

	###########################################################
	# extracting the reflection coefficient
	###########################################################
	yf = ypoints[-1]
	u = (n+1)*rh**(n+1)/(r(yf,n)**(n+2)*f(r(yf,n),n))
	ReB1 = u*sol[-1,2] - om*sol[-1,1]
	ImB1 = u*sol[-1,3] + om*sol[-1,0]
	ReB2 = u*sol[-1,2] + om*sol[-1,1]
	ImB2 = u*sol[-1,3] - om*sol[-1,0]

	refl = (ReB1**2 + ImB1**2)/(ReB2**2 + ImB2**2)
	return (1-refl)*degen(l,n,t)/gtot # returns the absorption coefficient, not greybody factor

#n = 0
#t = 1
#om = 4
#l = 2
#omeg = om*rh/(n+1)
#p0 = [A,0,0,-A*omeg]
#ypoints = ypo[n]
#sol = odeint(deriv, p0, ypoints, args=(om,l,n,t))
#plt.plot(ypoints, sol[:,0])
#plt.show()
#quit()

def cross(g, om, n): # turns an absorption coefficient into a greybody factor, units of rh**(n+2)
	#return pi*g/om**2 # the 4d result
	return 2**n*pi**((n+1)/2)*gamma((n+1)/2)*(n+1)*g/om**(n+2)

lmax = 30
lmin = 2
xpoints, xmin, xmax = 200, 0.01, 5
gs = []
gbs = []
xs = np.power(10, np.linspace(np.log10(xmin),np.log10(xmax),xpoints))
count = 0
outfile = "gb" + str(n) + str(s) + ".txt"
with open(outfile,'w') as myfile:
	for om in xs:
		count += 1
		g = 0
		for l in range(lmin,lmax):
			print(str(count)+"/"+str(xpoints)+' '+str(l) , end='\r', flush=True)
			for t in [0,1,2]:
				g += grey(om,l,n,t)
		gs.append(g)
		gbs.append(cross(g, om, n))
		myfile.write(str(g)+'\n')

plt.plot(xs, gbs)
plt.show()
