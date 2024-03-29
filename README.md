# MEF1D
A unidimensional finite element method (FEM) implemented in MATLAB to solve variational problems subject to boundary conditions.
## About
The finite element method (FEM) is a way to solve numerically differential equations by discretizing the function space where the solution lies. This particular implementation follows the approach presented in the chapter 5 of [La méthode des éléments finis: de la théorie à la pratique](https://dms.umontreal.ca/~owens/MAT6450/elfind.pdf) (A. Fortin & A. Garon). This code was made during an internship in applied maths.
## Usage
This code was made entirely for educational purposes, and in no means does it pretend to be fully optimized or efficient. The function *solve_mef1d()* was made to solve for *u* in unidimensional ODEs of order 2 with this particular form:
![Typical problem](https://github.com/paluneau/MEF1D/blob/master/problemetype.png)
where *p(x),q(x),R(x)* are parameters depending on the problem. Dirichlet conditions (*essential*) are imposed on the boundaries of the domain (*a,b*). Neumann conditions (*natural*) and mixed conditions (Dirichlet-Neumann) are also supported. **Linear** and **quadratic** shape functions (Lagrange) are available for the discretization. You can plot the resulting approximate solution with *plot_mef1d()*. Also, if you have computed the analytic solution to your problem and you just want to compare it with the numerical solution, *solve_mef1d()* can give you the L<sup>2</sup> error on your approximate solution, and *plot_mef1d()* can plot both on the same graph for a visual comparison.
## Validation
Code validity was ensured in 2 ways:
1. The Method of Manufactured Solution (MMS) was applied on a linear elasticity problem. Given an arbitrary analytical solution, the source term of the differential equation was computed by hand. By comparing the numerical solution (obtained with the hand-calculated source term) with the initial manufactured solution, we converged according to L<sup>2</sup> norm. It was also verified that linear solutions were exact while using order 1 and 2 elements, and that quadratic solutions were exact with order 2 elements.
2. Results from an example in the previously mentioned book were reproduced as expected.
## Future Functionalities
* An option to add "point loads" (Diracs) on specific nodes will be implemented eventually;
* Comments will be translated to english;
* Other example problems will be added;
* Maybe some refactoring of the code (dividing the code into smaller modules).

