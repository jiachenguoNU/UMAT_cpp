# UMAT_cpp guideline
User-defined material model in ABAQUS, coded in C++ (UMAT-C++)

J2 model with isotropic hardening. Power law is used for the hardening. Implemented using C++.
Eigen library is used for matrix and vector operation.

Since Eigen library is used for tensor computation, you need to compile the code into an obj file before using ABAQUS. 
To this end, you can use visual studio (2019 or later). When compile the obj file, make sure include the eigen library in your path. (You can download eigen library from their website.)

With the compiled obj file, you can run this umat with the command, `umat=j2-iso.obj` when you submit the job in the abaqus command window.
