################################################################################
#  Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

__precompile__()

"""
Local cubic interpolations in up to 4 dimensions. Interpolates
a given point using only local values and derivatives.
"""
module LocalInterpolations

import Base.diff
export Interpolation, interpolate, 
       diff, diff_x, diff_y, diff_z, diff_w, 
       Boundary

"""
Boundary condition for interpolation:
* `undef`:       Return NaN if out of bounds (default)
* `throw_error`: Throw an error if out of bounds
* `unsafe`:      Do not check indexes (dangerous)
"""
@enum Boundary undef throw_error unsafe #constant extrapolate periodic

#export Boundaries
for s in instances(Boundary)
    @eval export $(Symbol(s))
end 

"""
    Interpolation{T<:Number, B, N, M}

Interpolation coefficients type `T`, dimension `N` and size `M`, and boundary
condition B (see [`Boundary`](@ref)).
"""
struct Interpolation{T<:AbstractFloat, B, N, M}
    x0::NTuple{N, T}
    dx::NTuple{N, T}
    ks::Array{NTuple{M, T}, N} 
end

include("1D.jl")
include("2D.jl")
include("3D.jl")
include("4D.jl")

#in order not to wait 3min each time a 4D interpolation is used
precompile(cellcoeffs, ntuple(i->Array{Float32, 4}, 16))
precompile(cellcoeffs, ntuple(i->Array{Float64, 4}, 16))

for T in [Float32 Float64]
    for i in 0:1
     for j in 0:1
      for k in 0:1
       for l in 0:1
           precompile(poly, (NTuple{256, T}, ntuple(i->T, 4)..., 
                     dxdydzdw{i, j, k, l}))
        end
       end
      end
     end
end

#Functors
(ip::Interpolation)(x)          = interpolate(ip, x)
(ip::Interpolation)(x, y)       = interpolate(ip, x, y)
(ip::Interpolation)(x, y, z)    = interpolate(ip, x, y, z)
(ip::Interpolation)(x, y, z, w) = interpolate(ip, x, y, z, w)

################################################################################
#Private functions
################################################################################

nan(x::Rational)= 1//0
nan(x::Float16) = NaN16
nan(x::Float32) = NaN32
nan(x::Float64) = NaN64 
nan{T, N}(x::NTuple{N, T}) = ntuple(i -> nan(x[1]), N)

#return variable transformation and coefficient index
function xp_i{T<:AbstractFloat}(x::T, x0::T, dx::T)
    xr  = x - x0
    xt  = xr/dx

    ic  = floor(Int, xt) 
    xp  = xt - ic
    i   = ic + 1

    return xp, i
end

end
