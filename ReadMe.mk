#Elasto-Plastic Micromorphic Constitutive Model
email: nathanm@lanl.gov

Follows the derivations of Farhad Shahabi in his dissertation. Utilizes a 
general form of the elastic model derived from a quadratic Helmholtz free 
energy which would allow for anisotropy.

The plastic evolution equations are solved using a Newton-Raphson technique 
which can be modified to use a Homotopy Newton-Raphson solver if the problem 
is very stiff.
