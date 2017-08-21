#  Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, you can obtain one at http://mozilla.org/MPL/2.0/.

################################################################################
# 4D Interpolation
################################################################################

function Interpolation{T}(xs::Union{Array{T, 1}, StepRangeLen{T}}, 
                          ys::Union{Array{T, 1}, StepRangeLen{T}}, 
                          zs::Union{Array{T, 1}, StepRangeLen{T}}, 
                          ws::Union{Array{T, 1}, StepRangeLen{T}}, 
                          fs::Array{T, 4}, 
                          fs_x::Array{T, 4}, fs_y::Array{T, 4}, 
                          fs_z::Array{T, 4}, fs_w::Array{T, 4}, 
                          fs_xy::Array{T, 4}, fs_xz::Array{T, 4}, 
                          fs_xw::Array{T, 4}, 
                          fs_yz::Array{T, 4}, fs_yw::Array{T, 4}, 
                          fs_zw::Array{T, 4}, fs_xyz::Array{T, 4}, 
                          fs_xyw::Array{T, 4}, fs_xzw::Array{T, 4}, 
                          fs_yzw::Array{T, 4}, fs_xyzw::Array{T, 4}, 
                          b::Boundary=undef)

    x0  = (xs[1], ys[1], zs[1], ws[1])
    dx  = (xs[2] - xs[1], ys[2] - ys[1], zs[2] - zs[1], ws[2] - ws[1])
    sz  = size(fs)
    ks  = Array{NTuple{256, T}}(sz[1] - 1, sz[2] - 1, sz[3] - 1, sz[4] - 1)

    for i = 1:sz[1] - 1
        for j = 1:sz[2] - 1
            for k = 1:sz[3] - 1
                for l = 1:sz[4] - 1
                  ks[i, j, k, l] = 
                    cellcoeffs(
                    fs[i:i+1, j:j+1, k:k+1, l:l+1], 
                    fs_x[i:i+1, j:j+1, k:k+1, l:l+1]*dx[1],
                    fs_y[i:i+1, j:j+1, k:k+1, l:l+1]*dx[2], 
                    fs_z[i:i+1, j:j+1, k:k+1, l:l+1]*dx[3], 
                    fs_w[i:i+1, j:j+1, k:k+1, l:l+1]*dx[4], 
                    fs_xy[i:i+1, j:j+1, k:k+1, l:l+1]*dx[1]*dx[2],
                    fs_xz[i:i+1, j:j+1, k:k+1, l:l+1]*dx[1]*dx[3],
                    fs_xw[i:i+1, j:j+1, k:k+1, l:l+1]*dx[1]*dx[4],
                    fs_yz[i:i+1, j:j+1, k:k+1, l:l+1]*dx[2]*dx[3],
                    fs_yw[i:i+1, j:j+1, k:k+1, l:l+1]*dx[2]*dx[4],
                    fs_zw[i:i+1, j:j+1, k:k+1, l:l+1]*dx[3]*dx[4],
                    fs_xyz[i:i+1, j:j+1, k:k+1, l:l+1]*dx[1]*dx[2]*dx[3],
                    fs_xyw[i:i+1, j:j+1, k:k+1, l:l+1]*dx[1]*dx[2]*dx[4],
                    fs_xzw[i:i+1, j:j+1, k:k+1, l:l+1]*dx[1]*dx[3]*dx[4],
                    fs_yzw[i:i+1, j:j+1, k:k+1, l:l+1]*dx[2]*dx[3]*dx[4],
                    fs_xyzw[i:i+1, j:j+1, k:k+1, l:l+1]*dx[1]*dx[2]*dx[3]*dx[4])
                end
            end
        end
    end
    Interpolation{T, b, 4, 256}(x0, dx, ks)
end   


function Interpolation{T}(xs::Union{Array{T, 1}, StepRangeLen{T}}, 
                                ys::Union{Array{T, 1}, StepRangeLen{T}}, 
                                zs::Union{Array{T, 1}, StepRangeLen{T}}, 
                                ws::Union{Array{T, 1}, StepRangeLen{T}}, 
                                f::Function, f_x::Function,
                                f_y::Function, f_z::Function, 
                                f_w::Function, f_xy::Function, 
                                f_xz::Function, f_xw::Function, 
                                f_yz::Function, f_yw::Function, 
                                f_zw::Function, f_xyz::Function, 
                                f_xyw::Function, f_xzw::Function, 
                                f_yzw::Function, f_xyzw::Function, 
                                b::Boundary=undef)

    yns          = ys'
    
    zns          = Array{T}(1, 1, length(zs))
    zns[1, 1, :] = zs

    wns          = Array{T}(1, 1, 1, length(ws))
    wns[1, 1, 1, :] = ws

    Interpolation(xs, ys, zs, ws,  f.(xs, yns, zns, wns),     
                  f_x.(xs, yns, zns, wns), f_y.(xs, yns, zns, wns), 
                  f_z.(xs, yns, zns, wns), f_w.(xs, yns, zns, wns), 
                  f_xy.(xs, yns, zns, wns), f_xz.(xs, yns, zns, wns),
                  f_xw.(xs, yns, zns, wns), f_yz.(xs, yns, zns, wns),
                  f_yw.(xs, yns, zns, wns), f_zw.(xs, yns, zns, wns),
                  f_xyz.(xs, yns, zns, wns), f_xyw.(xs, yns, zns, wns),
                  f_xzw.(xs, yns, zns, wns), f_yzw.(xs, yns, zns, wns),
                  f_xyzw.(xs, yns, zns, wns), b)
end

function interpolate{T}(f::Function, ip::Interpolation{T, unsafe}, 
            x::T, y::T, z::T, w::T)
    xp, i = xp_i(x, ip.x0[1], ip.dx[1])
    yp, j = xp_i(y, ip.x0[2], ip.dx[2])
    zp, k = xp_i(z, ip.x0[3], ip.dx[3])
    wp, l = xp_i(w, ip.x0[4], ip.dx[4])

    @inbounds return f(ip.ks[i, j, k, l], xp, yp, zp, wp)
end

function interpolate{T}(f::Function, ip::Interpolation{T, throw_error}, 
                        x::T, y::T, z::T, w::T)
    xp, i = xp_i(x, ip.x0[1], ip.dx[1])
    yp, j = xp_i(y, ip.x0[2], ip.dx[2])
    zp, k = xp_i(z, ip.x0[3], ip.dx[3])
    wp, l = xp_i(w, ip.x0[4], ip.dx[4])

    return f(ip.ks[i, j, k, l], xp, yp, zp, wp)
end

function interpolate{T}(f::Function, ip::Interpolation{T, undef}, 
                        x::T, y::T, z::T, w::T)
    xp, i = xp_i(x, ip.x0[1], ip.dx[1])
    yp, j = xp_i(y, ip.x0[2], ip.dx[2])
    zp, k = xp_i(z, ip.x0[3], ip.dx[3])
    wp, l = xp_i(w, ip.x0[4], ip.dx[4])

    if checkbounds(Bool, ip.ks, i, j, k, l)
        @inbounds return f(ip.ks[i, j, k, l], xp, yp, zp, wp)
    else
        return nan(xp)
    end
end

interpolate(ip, x, y, z, w) = interpolate(poly, ip, x, y, z, w)

function diff(ip, x, y, z, w, edx, edy, edz, edw)
    poly_diff(ks, x, y, z, w) = poly(ks, x, y, z, w, 
                                dxdydzdw{edx, edy, edz, edw}())
    return interpolate(poly_diff, ip, x, y, z, w
           )/ip.dx[1]^edx/ip.dx[2]^edy/ip.dx[3]^edz/ip.dx[4]^edw
end

diff_x(ip, x, y, z, w)  = diff(ip, x, y, z, w, 1, 0, 0, 0)
diff_y(ip, x, y, z, w)  = diff(ip, x, y, z, w, 0, 1, 0, 0)
diff_z(ip, x, y, z, w)  = diff(ip, x, y, z, w, 0, 0, 1, 0)
diff_w(ip, x, y, z, w)  = diff(ip, x, y, z, w, 0, 0, 0, 1)

################################################################################
#Private functions
################################################################################

xterms  = [
    :(xt0=1; xt1=x; xt2=x*x; xt3=xt2*x),
    :(xt0=0; xt1=1; xt2=2*x; xt3=3*x*x),
    :(xt0=0; xt1=0; xt2=2;   xt3=6x),
    :(xt0=0; xt1=0; xt2=0;   xt3=6)
]

yterms  = [
    :(yt0=1; yt1=y; yt2=y*y; yt3=yt2*y),
    :(yt0=0; yt1=1; yt2=2*y; yt3=3*y*y),
    :(yt0=0; yt1=0; yt2=2;   yt3=6y),
    :(yt0=0; yt1=0; yt2=0;   yt3=6)
]

zterms  = [
    :(zt0=1; zt1=z; zt2=z*z; zt3=zt2*z),
    :(zt0=0; zt1=1; zt2=2*z; zt3=3*z*z),
    :(zt0=0; zt1=0; zt2=2;   zt3=6z),
    :(zt0=0; zt1=0; zt2=0;   zt3=6)
]

wterms  = [
    :(wt0=1; wt1=w; wt2=w*w; wt3=wt2*w),
    :(wt0=0; wt1=1; wt2=2*w; wt3=3*w*w),
    :(wt0=0; wt1=0; wt2=2;   wt3=6w),
    :(wt0=0; wt1=0; wt2=0;   wt3=6)
]

type dxdydzdw{M, N, O, P}
end

@generated function poly{T<:Number, M, N, O, P}(ks::NTuple{256 , T}, 
                        x::T, y::T, z::T, w::T, 
                        trait::dxdydzdw{M, N, O, P}=dxdydzdw{0, 0, 0, 0}())
    xt  = xterms[M+1]
    yt  = yterms[N+1]
    zt  = zterms[O+1]
    wt  = wterms[P+1]

    e   = :0
    for i = 0:3
        xi  = parse("xt$i")
        for j = 0:3
            yj  = parse("yt$j")
            for k = 0:3
                zk  = parse("zt$k")
                for l = 0:3
                    wl  = parse("wt$l")
                    e   = :($xi*$yj*$zk*$wl*ks[$i*64 + $j*16 + $k*4 + $l + 1] + 
                           $e)
                end
            end
        end
    end
    return :($xt; $yt; $zt; $wt; $e)
end

include("4Dcellcoeffs.jl")
