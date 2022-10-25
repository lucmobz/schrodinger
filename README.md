# Time Dependent Schrodinger Equation Solver

* Initial condition is a Gaussian wave packet
* Finite difference schemes in space (centered 2nd order) and time (Crank-Nicolson 2nd order)
* Uses Eigen library v3.4.0 for linear algebra and solvers
* Uses gnuplot for plotting
* Uses BiCGSTAB for solving the linear systems
* Units are such that h/2pi is 1 and x'=x/L, t'=t/(mL^2), V'=VmL^2 (L length, m mass, V potential)
