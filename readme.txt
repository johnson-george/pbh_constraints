tests.py is a script to deduce the format of the greybody tables that are required by blackhawk.

gb_tests.py is a first attempt at producing some greybody factors in the scalar case. In particular, it uses a different coordinate transformation (i.e. different coordinate y) for the ODE than that used in 4D.

gb_PARTICLE.py are four scripts that produce greybody factors in arbitrary spacetime dimension n, for scalars, fermions, gauge bosons, and gravitons.

format.py takes the output of gb_PARTICLE.py (a list of absorption coefficients) and converts it into a format that blackhawk can read. Note that it writes the same greybody factor for each choice of a = J/M, even though all computations involve a = 0.

photons.py does much of the heavy lifting. It takes the data from alldata.txt to determine the evolution of the black hole. It then computes the number of photons such a black hole would give rise to today.
