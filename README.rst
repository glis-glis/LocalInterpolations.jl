======================
LocalInterpolations.jl
======================

|travis-ci|

This package is a julia library for cubic, local interpolations in up to 4
dimensions. In order to calculate the coefficients of the cubic polynom, only
local values are used: The data itself and its first-order derivatives. This is
in contrast to splines, where the coefficients are not calculated using
derivatives, but non-local data, which can leed to over-smoothing the result.

.. |travis-ci| image:: https://api.travis-ci.org/AFueglistaler/LocalInterpolations.jl.svg?branch=master
    :target: https://travis-ci.org/AFueglistaler/LocalInterpolations.jl

Installation
============

To install the package, type

.. code:: 

    Pkg.clone("https://github.com/AFueglistaler/LocalInterpolations.jl.git")
    using LocalInterpolations

The second command precompiles the interpolations in 4D and will take a while
(5-10 minutes).

Usage
=====
1D
__

.. code:: julia

    #Define function and derivative
    f(x) = log(x^2 + 1)
    f_x(x) = 2x/(x^2 + 1)
    #Define interpolation range
    xs = 1:0.1:3
    #Create interpolation
    ip = Interpolation(xs, f, f_x)
    #interpolate at 1.33 and compare with exact value
    ip(1.33) - f(1.33)
    #interpolate the derivative at 1.33 and compare with exact value
    diff_x(ip, 1.33) - f_x(1.33)
    #second order derivative
    diff(ip, 1.33, 2)
    #third order derivative
    diff(ip, 1.33, 3)
    #There is no derivative for cubic polynoms higher than 3!

There are different behaviours outside the interpolation range: 

+ ``undef`` returns nan outside of range (default)
+ ``throw_error`` throws an error
+ ``unsave`` does nothing (fast but dangerous)
   
You can define the boundary behaviour as an additional argument:
   
.. code:: julia
      
    julia> ip = Interpolation(xs, f, f_x, undef)
    julia> ip(10.)
        NaN
    julia> ip = Interpolation(xs, f, f_x, throw_error)
    julia> ip(10.)
        ERROR: BoundsError: attempt to access 20-element Array{NTuple{4,Float64},1} at index [90]
    julia> ip = Interpolation(xs, f, f_x, unsafe)
    julia>#Do not do this! ip(10.)
   
Higher dimensions
-----------------
Using LocalInterpolations with higher dimensions is the same as for 1D, but you
need many more derivatives. For example, in 3 dimensions, you will need
:math:`f, f_x =  df/dx, f_y, f_z, f_xy = d²f/dx/dy, f_xz, f_yz` and
:math:`f_xyz`. 

.. code:: julia

    #Define function and derivative (all the same)
    exp4d(x, y, z, w) = exp(x + y + z + w)
    #Define interpolation range
    xs = 1:0.1:2
    #Create interpolation
    ip = Interpolation(ntuple(i->xs, 4)..., ntuple(i->exp4d, 16)...)
    #interpolate at (1.13, 1.24, 1.35, 1.46) and compare with exact value
    ip(1.13, 1.24, 1.35, 1.46) - exp4d(1.13, 1.24, 1.35, 1.46)
    #interpolate the derivative in z at the same point and compare
    diff_z(ip, 1.13, 1.24, 1.35, 1.46) - exp4d(1.13, 1.24, 1.35, 1.46)
