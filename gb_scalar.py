import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
import sys

n = int(sys.argv[1])
s = 0
rh = 1

def f(r,n): 
	return 1 - (rh/r)**(n+1)

def rr(y,n):
	if y < -10**(-10):
		out = rh*(1-np.exp(y))**(-1/(n+1))
	else:
		out = rh*(-y)**(-1/(n+1))
	return out
r = np.vectorize(rr)

def deriv(p,y,om,l,n):
	ReR, ImR, ReS, ImS = p
	k = (n+1)*rh**(n+1)
	h = f(r(y,n),n)
	s = r(y,n)
	a = h**2*s**(2*n+2)*n*(n/2+1)/(2*k**2)
	b = h*s**(n+1)*n/(2*k)
	c = s**(2*n+4)*om**2/k**2
	d = -s**(2*n+2)*h*l*(l+1)/k**2
	x = -a-b-c-d
	return [ReS, ImS, x*ReR, x*ImR]

###########################################################
# initialising the solver 
###########################################################
steps = 10**6
A = 1
ypo = []
for m in range(7):
	points = np.linspace(np.log10(5*(m+1)), -3*(m+1), steps+1)
	ypo.append(-np.power(10, points))

def grey(om,l,n):
	###########################################################
	# solving the ODE
	###########################################################
	omeg = om*rh/(n+1)
	p0 = [A,0,0,-A*omeg]
	ypoints = ypo[n]

	#yspan = (ypoints[0], ypoints[-1])
	#soln = solve_ivp(lambda y,p: deriv(p,y,om,l,n), yspan, p0, t_eval=ypoints)
	#sol = soln.y
	sol = odeint(deriv, p0, ypoints, args=(om,l,n))

	###########################################################
	# extracting the reflection coefficient
	###########################################################
	yf = ypoints[-1]
	u = (n+1)*rh**(n+1)/(r(yf,n)**(n+1)*f(r(yf,n),n))
	m = n/2+1
	ReB1e = u*sol[-1,2] + m*sol[-1,0] - om*sol[-1,1]*r(yf,n)
	ImB1e = u*sol[-1,3] + m*sol[-1,1] + om*sol[-1,0]*r(yf,n)
	ReB2e = u*sol[-1,2] + m*sol[-1,0] + om*sol[-1,1]*r(yf,n)
	ImB2e = u*sol[-1,3] + m*sol[-1,1] - om*sol[-1,0]*r(yf,n)

	c = np.cos(om*r(yf,n))
	s = np.sin(om*r(yf,n))
	ReB1 = c*ReB1e + s*ImB1e
	ImB1 = c*ImB1e - s*ReB1e
	ReB2 = c*ReB2e - s*ImB2e
	ImB2 = c*ImB2e + s*ReB2e

	refl = (ReB1**2 + ImB1**2)/(ReB2**2 + ImB2**2)
	GB = np.pi*(1-refl)*(2*l+1)/om**2
	#return GB/(rh**2*np.pi) 
	return (1-refl)*(2*l+1) # returns the absorption coefficient, not greybody factor

#n = 0
#om = 0.36
#l = 0
#omeg = om*rh**2
#p0 = [A,0,0,-A*omeg]
#ypoints = ypo[n]
#sol = odeint(deriv, p0, ypoints, args=(om,l,n))
#plt.plot(ypoints, sol[:,0])
#plt.show()
#quit()

lmax = 30
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
		for l in range(0,lmax):
			print(str(count)+"/"+str(xpoints)+' '+str(l) , end='\r', flush=True)
			g += grey(om,l,n)
		gs.append(g)
		gbs.append(g/(rh**2*om**2))
		myfile.write(str(g)+'\n')
		
#gbs = []
#with open(outfile, 'r') as myfile:
#	for line, om in zip(myfile, xs):
#		gb = float(line[0:-1])
#		gbs.append(gbb/(rh**2*om**2))

plt.plot(xs, gbs)
plt.show()
