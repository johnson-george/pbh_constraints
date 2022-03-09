import math
import numpy as np
from scipy.integrate import solve_ivp, quad
import matplotlib.pyplot as plt
from gammaraydata import data_test, data

###########################################################
# create a table of emission rates from alldata.txt
###########################################################

dPdt = [] # number of photons produced per time per volume as a function of (M, E)
dNdt = [] # number of particles produced per time per volume as a function of (M, E)

mmin, mmax, mtot = 9, 19, 1001
masses = np.power(10, np.linspace(mmin, mmax, mtot))
emin, emax, etot = 1e-10, 1e+10, 2000 # the bounds on energy in the primary spectra
emin2, emax2, etot2 = 1e-6, 1e+5, 500 # the bounds on energy in the secondary spectra
energies = np.exp(np.linspace(np.log(emin), np.log(emax), etot))
energies2 = np.exp(np.linspace(np.log(emin2), np.log(emax2), etot2))

print("reading data files...")

with open("alldata.txt", "r") as myfile:
	dat = myfile.readlines()
	dat = [line.split()[1:] for line in dat]
	for n in range(mtot):
		start = etot*n + n + 1
		mdata = dat[start:start + etot]
		N = [sum([float(partic) for partic in entry]) for entry in mdata]
		dNdt.append(np.array(N))

with open("photdata.txt", "r") as myfile2:
	dat = myfile2.readlines()
	dat = [line.split()[1:] for line in dat]
	for n in range(mtot):
		start = etot2*n + n + 1
		mdata = dat[start:start + etot]
		P = [float(entry[0]) for entry in mdata]
		dPdt.append(np.array(P))

print("reading: complete")

###########################################################
# first we calculate M(t), the mass as a function of time
###########################################################

def spec(M): # this function linearly interpolates the above data to arbitrary M
	step = (mmax - mmin)/(mtot - 1)	
	expon = np.log10(M)
	d = expon % step
	r = expon - d
	i = int(round((r - mmin)/step))
	return (1 - d/step)*dNdt[i] + (d/step)*dNdt[i+1]

def totalrate(spec): # this function essentially integrates over energy
	dEdt = np.multiply(spec, energies)
	E = 0
	for j in range(etot - 1):
		E += 0.5*(dEdt[j+1] + dEdt[j])*(energies[j+1] - energies[j])	
	return E

def power(t, M, start): # returns dM/dt
	M = M[0]
	conv = 1.782662e-24
	if M < start/20 or M < 10**mmin:
		X = 0
	else:
		X = -conv*totalrate(spec(M))	
	return X

print("calculating evolution of black hole...")

lifetimes = [4e+2,5e+5,5e+8,7e+11,8e+14,1.2e+18,4e+21,6e+24,6e+27,6e+30] # lifetime guesses
crittimes = [] # the lifetime of each order of magnitude
M0 = 10**mmax
finaltimes = [] # this contains 10 arrays of times, all starting at zero
finalsol = [] # this contains 10 arrays of masses, with finalsol[n][-1] = finalsol[n+1][0]
for m in reversed(range(mmax-mmin)):
	tf = 2*lifetimes[m]
	teval = np.linspace(0, tf, 10**6)
	def wrapper(t, M):
    		return power(t, M, M0)
	M = solve_ivp(wrapper, (0, tf), np.array([M0]), t_eval=teval)
	Ms = [x for x in M.y[0] if x > 10**(mmin+m)]
	ts = M.t[0:len(Ms)]
	M0 = Ms[-1]
	finaltimes.append(ts)
	finalsol.append(Ms)
	crittimes.append(ts[-1])

print("evolution: complete")

###########################################################
# next we calculate I(E), intensity as a function of energy
###########################################################

t0 = 4.35*10**17 # time today
tmin = 1.20*10**13 #380,000 years
deltat = t0 - tmin

lifetimes = [sum(crittimes[k:]) for k in range(10)] # redefining the exact lifetime
ll = len([x for x in lifetimes if x > deltat]) # number of black holes older than the universe
remainder = deltat - lifetimes[ll]
indexc = next(i for i,time in enumerate(finaltimes[ll-1]) if time>crittimes[ll-1]-remainder) 
tc = finaltimes[ll-1][indexc-1]
Mc = finalsol[ll-1][indexc-1] # the initial black hole mass which evaporates today
print("critical mass = ", Mc)

def z(t): # redshift --- note discrepancy with Silk paper!
	return (t0/t)**(2/3) - 1
zmax = z(tmin)

def mcont(t, a, b): # creating the mass interpolation, evolving Mi from time 0 to time t
	j = next(i for i,s in enumerate(b) if s>t) # the first entry where time is more than t
	return a[j-1] + (t - b[j-1])*(a[j] - a[j-1])/(b[j] - b[j-1])

def specp(M): # another interpolater, but just for photons
	step = (mmax - mmin)/(mtot - 1)	
	expon = np.log10(M)
	d = expon % step
	r = expon - d
	i = int(round((r - mmin)/step))
	return (1 - d/step)*dPdt[i] + (d/step)*dPdt[i+1]

def phot(t, E, a, b): # the integrand to be integrated over time
	M = mcont(t, a, b)
	S = specp(M)
	Ethen = E*(1 + z(t))
	j = next(i for i,f in enumerate(energies2) if f>Ethen)
	intp = S[j-1] + (Ethen - energies2[j-1])*(S[j] - S[j-1])/(energies2[j] - energies2[j-1])
	return intp*(1 + z(t)) 

c0 = 29979245800 # speed of light in cm per s
def I(E, a, b, tmax):
	N = quad(phot, tmin, tmax, args=(E, a, b,), limit=50)[0] # number of photons per GeV
	lum = N*E
	return c0*lum/(4*math.pi)

###########################################################
# comparing the intensity with constraints, for certain Mi
###########################################################

rhodm = 2.25e-30 # dark matter density in g cm-3

critmasses = [m[0] for m in finalsol]
def part_evolution(Mi): # this function takes the initial mass and gives the evolution from then
	k = next(i for i,m in enumerate(critmasses) if m<Mi) - 1
	kk = next(i for i,m in enumerate(finalsol[k]) if m<Mi) - 1 # mass just bigger than Mi
	apre = finalsol[k][kk:]
	a = apre + [item for sublist in finalsol[k+1:] for item in sublist]
	b = [time - finaltimes[k][kk] for time in finaltimes[k][kk:]]
	for p in range(9-k):
		bf = b[-1]
		bpre = [time + bf for time in finaltimes[k+p+1]]
		b.extend(bpre)
	tfin = b[-1] # the lifetime of the black hole
	return (a, b, tfin)

f = []
massrange = np.linspace(13, 17, 81)
for n in massrange:
	Mi = 10**n
	print(n)
	(a, b, tfin) = part_evolution(Mi)
	if tfin < tmin:
		f.append(10**10)
		continue
	tmax = min(t0, tfin) # stop integrating either today or when the black hole is gone
	cons = []
	checkable_energies = [E for E in energies2 if E < energies2[-1]/zmax]
	
	for E in checkable_energies:
		F = I(E, a, b, tmax)
		if F <= 0: 
			continue
		comp = np.log10(F)
		expr = np.log10(data(E))
		logrho = expr - comp
		cons.append(logrho)

	fm = (10**min(cons))*Mi/rhodm
	f.append(fm)

print(f)
with open("constraints.txt", "w") as myfile:
	for entry in f:
		myfile.write(str(entry))
		myfile.write("\n")
plt.plot(massrange, np.log10(f))
plt.show()

omegac = 0.264
h0 = 67.4
