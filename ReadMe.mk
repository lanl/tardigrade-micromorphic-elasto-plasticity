#Elasto-Plastic Micromorphic Constitutive Model
email: nathanm@lanl.gov

C20048 Tardigrade

 2021. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.

Follows the derivations of Farhad Shahabi in his dissertation. Utilizes a 
general form of the elastic model derived from a quadratic Helmholtz free 
energy which would allow for anisotropy.

The plastic evolution equations are solved using a Newton-Raphson technique 
which can be modified to use a Homotopy Newton-Raphson solver if the problem 
is very stiff.
