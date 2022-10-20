# Top88 MMA

Implementation of the Top88 Matlab code [1] with small changes to run a simplified version of the MMA method [2] for minimum compliance topology optimization. 

### Simplified MMA

Due to the derivative of the compliance always being negative (dc/dx < 0) the p term of the MMA function f(x) is always zero. Thus the MMA function f(x) can be simplified into the following expression:

f(x) = r + Σ(q/(x-L)

Using the simplified form of the MMA function the Lagrangian function is described using the following expression: 

L(x,λ) = f(x) + λ (Σx - volfrac)

To compute the new set of material densities x the derivative of the Lagrangian function with respect to x must be set equal to zero following the expression:

dL/dx = q/(x-L)^2 + λdv = 0 => x = sqrt(q/λdv) + L

Thus the final expression for computing the new set of material densities x is following the expression:

xnew = sqrt(q/λdv) + L

Where L is one of the two moving asymptotes updated using the expressions defined in [2] and term q is taken from the MMA function definition:

q = -(x-L)^2 dC/dx

[1] E. Andreassen, A. clausen, M. Schevenels, B.S. Lazarov, O. Sigmund. Efficient topology optimization in MATLAB using 88 lines of code, _structural multidisciplinary optimization_, 2011. 

[2] Kr. Svanberg. The method of moving asymptotes - a new method for structural optimization, _International Journal for numerical methods in Engineering_, 1987.
