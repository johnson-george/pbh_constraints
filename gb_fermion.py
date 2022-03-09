import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
import sys

n = int(sys.argv[1])
s = 0.5
rh = 1

def f(r,n): 
	return 1 - (rh/r)**(n+1)

def deriv(p,r,om,l,n):
	ReR, ImR, ReS, ImS = p
	h = f(r,n)
	der = 1/r + (n+1)*(1-h)/(2*h*r) 
	img = om*(1 - (n+1)*(1-h)/(2*h))/(r*h)
	fre = om**2/h**2
	ang = (l*(l+1) + 0.25)/(h*r**2)
	return [ReS, ImS, -der*ReS + (ang-fre)*ReR + img*ImR, -der*ImS + (ang-fre)*ImR - img*ReR]

###########################################################
# initialising the solver
###########################################################
steps = 10**6
A = 1
eps = 10**(-5)
h0 = f(1+eps,n)

def grey(om,l,n):
	###########################################################
	# solving the ODE
	###########################################################
	rpoints = np.linspace(1+eps, int(10**4/om), steps)
	omeg = om/h0
	p0 = [A,0,0,-A*omeg]
	sol = odeint(deriv, p0, rpoints, args=(om,l,n))

	###########################################################
	# extracting the reflection coefficient
	###########################################################
	rf = rpoints[-1]
	c = np.cos(om*rf)
	s = np.sin(om*rf)
	ReP = sol[-1,0]
	ImP = sol[-1,1]

	ReB = ReP*c - ImP*s
	ImB = ReP*s + ImP*c
	absorp = A**2/(ReB**2 + ImB**2)
	return absorp*(2*l+1) # returns the absorption coefficient, not greybody factor

#n = 0
#om = 0.36
#j = 0.5
#omeg = om/h0
#p0 = [A,0,0,-A*omeg]
#sol = odeint(deriv, p0, rpoints, args=(om,j,n))
#plt.plot(rpoints, sol[:,0])
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
			j = l + 0.5
			print(str(count)+"/"+str(xpoints)+' '+str(l) , end='\r', flush=True)
			g += grey(om,j,n)
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
