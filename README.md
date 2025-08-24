Date: 8/23/2025

Instructions for using the 3d-RAIL source code.
Language: Originally written in MATLAB, version R2023b

If issues or questions arise, please email Joseph Nakao (Swarthmore College) at jnakao1@swarthmore.edu.

If this source code is used, please cite the original papers of the major algorithms utilized in the 3d-RAIL scheme:

* "A low-rank, high-order implicit-explicit integrator for three-dimensional convection-diffusion equations" by J. Nakao, G. Ceruti, and L. Einkemmer, (2025).
* "Numerical solution of a class of third order tensor linear equations" by V. Simoncini, (2020).
* "A Local Macroscopic Conservative (LoMaC) low rank tensor method for the Vlasov dynamics" by W. Guo and J.-M. Qiu, (2024).

### Step 1.

You must [download the Tensorlab package](https://tensorlab.net/). As stated in the Tensorlab documentation, please cite:

* "Tensorlab 3.0--numerical optimization strategies for large-scale constrained and coupled matrix/tensor factorization" by N. Vervliet, O. Debals, and L. De Lathauwer, (2016).
* "Tensorlab 3.0" by N. Vervliet et al., (2016).

After downloading the Tensorlab library, add it to the path.

### Step 2.

We have two separate folders: one for linear advection-diffusion equations, and another for the nonlinear viscous Burgers' equation. Both contain mostly the same files with only minor differences, e.g., for computing the flux terms. Regardless of which equation you solve, the following parameters must be changed in the main.m file:

* Final time, Tf
* Truncation tolerance, tol
* Mesh size, Nx, Ny, Nz
* CFL constant. If doing a single run, let Lambdavals be the desired constant. If testing the order of accuracy (as determined by compute_error=0,1), let Lambdavals be the vector of CFL constants.
* Scale for the LoMaC weight function, s
* Whether you want to use IMEX111, IMEX222, or IMEX443
* The numerical test from the 3d-RAIL paper to run, testnumber (see test_parameters.m)

### Step 3.

You must go into the IMEX111.m, IMEX222.m, and IMEX443.m files and manually (un)comment the truncation procedure at the end of each Runge-Kutta stage that you would like to use. There are four options: nonconservative, mass conservative, mass and momentum conservative, and mass and momentum and energy conservative. Be sure to change it at each place the truncation is performed, including in IMEX111 since it is used in the prediction step.

We've implemented a conditional loop in the main script main.m that pauses the simulation if the multilinear rank is getting too large. It is located in the time-stepping for-loop. Change it or comment it out, as desired.

We also compute the number density, bulk velocity, and energy (and L1 decay and relative entropy in the DFP equation) after each time-step. These are also located in the time-stepping for-loop. If any of these quantities are undesired, then comment them out to eliminate unnecessary computation. These are only intended for full verification of structure preservation.

If one wants to apply the 3d-RAIL algorithm to another similar PDE, then go into test_parameters.m and add another test. {A,B,C} construct the flow field, and {P1,P2,P3,Q} construct the source term. You also need to define the domain, CFLconstraints, and diffusion coefficients. In the current framework, we set parameters u_exact=0 and computer_error=0 if only a single run is being run and the error does not need to be computed. This can easily be changed if needed.
