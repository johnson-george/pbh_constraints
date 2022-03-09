import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

rh = 1
m = 0

def f(r,n): 
	return 1 - (rh/r)**(n+1)

def r(y,n):
	d4 = rh/(1-np.exp(rh*y))
	X = np.exp(2*rh*y)
	d5 = rh*(1+X)/(1-X)
	return n*d5 + (1-n)*d4

# presently takes no account of mass
def deriv(p,y,om,l,n):
	ReR, ImR, ReS, ImS = p
	ang = l*(l+1)/r(y,n)**2
	x = -(r(y,n)**4)*(om**2 - ang*f(r(y,n),n))
	return [ReS, ImS, x*ReR, x*ImR]

###########################################################
# initialising the solver
###########################################################
steps = 10**6
A = 1
ypoints = np.linspace(-3, -10**(-3), steps+1)

def grey(om,l,n):
	###########################################################
	# solving the ODE
	###########################################################
	omeg = om*rh**2
	p0 = [A,0,0,-A*omeg]
	sol = odeint(deriv, p0, ypoints, args=(om,l,n))

	###########################################################
	# extracting the reflection coefficient
	###########################################################
	yf = ypoints[-1]
	h = 1/(r(yf,n)*f(r(yf,n),n))
	ReB1e = h*sol[-1,2] + sol[-1,0] - om*sol[-1,1]*r(yf,n)
	ImB1e = h*sol[-1,3] + sol[-1,1] + om*sol[-1,0]*r(yf,n)
	ReB2e = h*sol[-1,2] + sol[-1,0] + om*sol[-1,1]*r(yf,n)
	ImB2e = h*sol[-1,3] + sol[-1,1] - om*sol[-1,0]*r(yf,n)

	c = np.cos(om*r(yf,n))
	s = np.sin(om*r(yf,n))
	ReB1 = c*ReB1e + s*ImB1e
	ImB1 = c*ImB1e - s*ReB1e
	ReB2 = c*ReB2e - s*ImB2e
	ImB2 = c*ImB2e + s*ReB2e

	refl = (ReB1**2 + ImB1**2)/(ReB2**2 + ImB2**2)
	GB = np.pi*(1-refl)*(2*l+1)/om**2
	return GB/(rh**2*np.pi)

#n = 0
#om = 0.36
#l = 0
#omeg = om*rh**2
#p0 = [A,0,0,-A*omeg]
#sol = odeint(deriv, p0, ypoints, args=(om,l,n))
#plt.plot(ypoints, sol[:,0])
#plt.show()
#quit()

oms = []
g4s = []
g5s = []
count = 0
lmax = 1
with open('trials.txt','w') as myfile:
	for om in np.linspace(0.01,5,100):
		count += 1
		print(str(count), end='\r', flush=True)
		g4 = 0
		g5 = 0		
		for l in range(0,lmax):
			g4 += grey(om,l,0)
			g5 += 1 #grey(om,l,1)
		oms.append(om)
		g4s.append(g4)
		g5s.append(g5)
		myfile.write(str(om)+" "+str(g4)+"\n")

plt.plot(oms, g4s)
#plt.plot(oms, g5s)
plt.show()
