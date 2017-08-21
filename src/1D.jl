#  Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

################################################################################
# 1D Interpolation
################################################################################

"""
    Interpolation(xs, fs, fs_x, b=undef)                            #1D

    Interpolation(xs, ys, fs, fs_x, fs_y, fs_xy, b=undef)           #2D

    Interpolation(xs, ys, zs, fs, fs_x, fs_y, fs_z, 
                  fs_xy, fs_xz, fs_yz, fs_xyz, b=undef)             #3D

    Interpolation(xs, ys, zs, ws, fs, fs_x, fs_y, fs_z, fs_w,
                  fs_xy, fs_xz, fs_xw, fs_yz, fs_yw, fs_zw, 
                  fs_xyz, fs_xyw, fs_xzw, fs_yzw, fs_xyzw, b=undef) #4D

Create interpolation in range[s] `xs`[, `ys`, `zs`, `ws`] for values `fs` =
f(`xs`[, `ys`, `zs`, `ws`]) and derivative values `fs_x` = df/dx(`xs`[, `ys`,
`zs`, `ws`])[, ..., `fs_xyzw`=d⁴f/(dx dy dz dw)(`xs`, `ys`, `zs`, `ws`)] and
boundary condition b.  
"""
function Interpolation{T}(xs::Union{Array{T, 1}, StepRangeLen{T}}, 
                                fs::Array{T, 1}, fs_x::Array{T, 1}, 
                                b::Boundary=undef)
    x0  = xs[1]
    dx  = xs[2] - xs[1]
    l   = length(fs) - 1 
    ks  = Array{NTuple{4, T}}(l)
    for i = 1:l
        ks[i] = cellcoeffs(fs[i:i+1], fs_x[i:i+1]*dx)
    end
    Interpolation{T, b, 1, 4}(tuple(x0), tuple(dx), ks)
end   

"""
    Interpolation(xs, f, f_x, b=undef)                          #1D

    Interpolation(xs, ys, f, f_x, f_y, f_xy, b=undef)           #2D

    Interpolation(xs, ys, zs, f, f_x, f_y, f_z, 
                  f_xy, f_xz, f_yz, f_xyz, b=undef)             #3D

    Interpolation(xs, ys, zs, ws, f, f_x, f_y, f_z, f_w,
                  f_xy, f_xz, f_xw, f_yz, f_yw, f_zw, 
                  f_xyz, f_xyw, f_xzw, f_yzw, f_xyzw, b=undef)  #4D

Create interpolation in range[s] `xs`[, `ys`, `zs`, `ws`] for function[s]
`f`(x[, y, z, w]) and derivative(s) `f_x` = df/dx(x[, y, z, w])[, ...,
`f_xyzw`=d⁴f/(dx dy dz dw)(x, y, z, w)] and boundary condition b.  
"""
function Interpolation{T}(xs::Union{Array{T, 1}, StepRangeLen{T}}, 
                                f::Function, f_x::Function, 
                                b::Boundary=undef)
    Interpolation(xs, f.(xs), f_x.(xs), b)
end   

function interpolate{T}(f::Function, ip::Interpolation{T, unsafe}, x::T)
    xp, i = xp_i(x, ip.x0[1], ip.dx[1])
    @inbounds return f(ip.ks[i], xp)
end

function interpolate{T}(f::Function, ip::Interpolation{T, throw_error}, x::T)
    xp, i = xp_i(x, ip.x0[1], ip.dx[1])
    return f(ip.ks[i], xp)
end

function interpolate{T}(f::Function, ip::Interpolation{T, undef}, x::T)
    xp, i = xp_i(x, ip.x0[1], ip.dx[1])

    if checkbounds(Bool, ip.ks, i)
        @inbounds return f(ip.ks[i], xp)
    else
        return nan(xp)
    end
end

"""
    interpolate(ip, x[, y, z, w])

interpolate with interpolation `ip` at point (`x`[, `y`, `z`, `w`]).
"""    
interpolate(ip, x)  = interpolate(poly, ip, x)

"""
    diff(ip, x[, y, z, w], edx[, edy, edz, edw])

interpolate the (`edx`[, `edy`, `edz`, `edw`])-order derivative with
interpolation `ip` at point (`x`[, `y`, `z`, `w`]). For example, interpolate in
2D the second order derivative in x and first order in y:

    julia> diff(ip, x0, y0, 2, 1)
"""    
function diff(ip, x, edx)
    poly_diff(ks, x)   = poly(ks, x, dx{edx}())
    return interpolate(poly_diff, ip, x)/ip.dx[1]^edx
end

"""
    diff_x(ip, x[, y, z, w])

First order derivative in x with interpolation `ip` at point (`x`[, `y`, `z`,
`w`]). 
""" 
diff_x(ip, x)   = diff(ip, x, 1)

"""
    diff_y(ip, x, y[, z, w])

First order derivative in y with interpolation `ip` at point (`x`, `y`[, `z`,
`w`]). 
""" 
diff_y

"""
    diff_z(ip, x, y, z[, w])

First order derivative in z with interpolation `ip` at point (`x`, `y`, `z`[,
`w`]). 
""" 
diff_z

"""
    diff_w(ip, x, y, z, w)

First order derivative in w with interpolation `ip` at point (`x`, `y`, `z`,
`w`). 
""" 
diff_w

################################################################################
#Private functions
################################################################################

xterms  = [
    :(xt0=1; xt1=x; xt2=x*x; xt3=xt2*x),
    :(xt0=0; xt1=1; xt2=2*x; xt3=3*x*x),
    :(xt0=0; xt1=0; xt2=2;   xt3=6x),
    :(xt0=0; xt1=0; xt2=0;   xt3=6)
]

type dx{N}
end

@generated function poly{T<:Number, N}(ks::NTuple{4 , T}, x::T, 
           trait::dx{N}=dx{0}())
    xt  = xterms[N+1]

    e   = :0
    for i = 0:3
        xi  = parse("xt$i")
        e   = :($xi*ks[$i + 1] + $e)
    end
    return :($xt; $e)
end

#coefficients for one cell
cellcoeffs{T<:Number}(fs::Array{T, 1}, fs_x::Array{T, 1}) = (
   fs[1],
   fs_x[1],
   -3fs[1] + 3fs[2] - 2fs_x[1] - fs_x[2],
   2fs[1] - 2fs[2] + fs_x[1] + fs_x[2]
)
