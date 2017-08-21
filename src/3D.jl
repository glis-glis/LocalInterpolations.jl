#  Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, you can obtain one at http://mozilla.org/MPL/2.0/.

################################################################################
# 3D Interpolation
################################################################################

function Interpolation{T}(xs::Union{Array{T, 1}, StepRangeLen{T}}, 
                                ys::Union{Array{T, 1}, StepRangeLen{T}}, 
                                zs::Union{Array{T, 1}, StepRangeLen{T}}, 
                                fs::Array{T, 3}, fs_x::Array{T, 3}, 
                                fs_y::Array{T, 3}, fs_z::Array{T, 3}, 
                                fs_xy::Array{T, 3}, fs_xz::Array{T, 3}, 
                                fs_yz::Array{T, 3}, fs_xyz::Array{T, 3}, 
                                b::Boundary=undef)
    x0  = (xs[1], ys[1], zs[1])
    dx  = (xs[2] - xs[1], ys[2] - ys[1], zs[2] - zs[1])
    sz  = size(fs)
    ks  = Array{NTuple{64, T}}(sz[1] - 1, sz[2] - 1, sz[3] - 1)
    for i = 1:sz[1] - 1
        for j = 1:sz[2] - 1
            for k = 1:sz[3] - 1
                ks[i, j, k] = cellcoeffs(fs[i:i+1, j:j+1, k:k+1], 
                            fs_x[i:i+1, j:j+1, k:k+1]*dx[1],
                            fs_y[i:i+1, j:j+1, k:k+1]*dx[2], 
                            fs_z[i:i+1, j:j+1, k:k+1]*dx[3], 
                            fs_xy[i:i+1, j:j+1, k:k+1]*dx[1]*dx[2],
                            fs_xz[i:i+1, j:j+1, k:k+1]*dx[1]*dx[3],
                            fs_yz[i:i+1, j:j+1, k:k+1]*dx[2]*dx[3],
                            fs_xyz[i:i+1, j:j+1, k:k+1]*dx[1]*dx[2]*dx[3])
            end
        end
    end
    Interpolation{T, b, 3, 64}(x0, dx, ks)
end   

function Interpolation{T}(xs::Union{Array{T, 1}, StepRangeLen{T}}, 
                                ys::Union{Array{T, 1}, StepRangeLen{T}}, 
                                zs::Union{Array{T, 1}, StepRangeLen{T}}, 
                                f::Function, f_x::Function,
                                f_y::Function, f_z::Function, 
                                f_xy::Function, f_xz::Function, 
                                f_yz::Function, f_xyz::Function, 
                                b::Boundary=undef)

    yns          = ys'
    
    zns          = Array{T}(1, 1, length(zs))
    zns[1, 1, :] = zs

    Interpolation(xs, ys, zs,  f.(xs, yns, zns),     
                  f_x.(xs, yns, zns), f_y.(xs, yns, zns), 
                  f_z.(xs, yns, zns), f_xy.(xs, yns, zns), 
                  f_xz.(xs, yns, zns),f_yz.(xs, yns, zns),
                  f_xyz.(xs, yns, zns), b)
end



function interpolate{T}(f::Function, ip::Interpolation{T, unsafe}, 
            x::T, y::T, z::T)
    #using varargs leads to performance issues
    xp, i = xp_i(x, ip.x0[1], ip.dx[1])
    yp, j = xp_i(y, ip.x0[2], ip.dx[2])
    zp, k = xp_i(z, ip.x0[3], ip.dx[3])

    @inbounds return f(ip.ks[i, j, k], xp, yp, zp)
end

function interpolate{T}(f::Function, ip::Interpolation{T, throw_error}, 
                        x::T, y::T, z::T)
    xp, i = xp_i(x, ip.x0[1], ip.dx[1])
    yp, j = xp_i(y, ip.x0[2], ip.dx[2])
    zp, k = xp_i(z, ip.x0[3], ip.dx[3])

    return f(ip.ks[i, j, k], xp, yp, zp)
end

function interpolate{T}(f::Function, ip::Interpolation{T, undef}, 
                        x::T, y::T, z::T)
    xp, i = xp_i(x, ip.x0[1], ip.dx[1])
    yp, j = xp_i(y, ip.x0[2], ip.dx[2])
    zp, k = xp_i(z, ip.x0[3], ip.dx[3])

    if checkbounds(Bool, ip.ks, i, j, k)
        @inbounds return f(ip.ks[i, j, k], xp, yp, zp)
    else
        return nan(xp)
    end
end

interpolate(ip, x, y, z) = interpolate(poly, ip, x, y, z)

function diff(ip, x, y, z, edx, edy, edz)
    poly_diff(ks, x, y, z) = poly(ks, x, y, z, dxdydz{edx, edy, edz}())
    return interpolate(poly_diff, ip, x, y, z
           )/ip.dx[1]^edx/ip.dx[2]^edy/ip.dx[3]^edz
end

diff_x(ip, x, y, z)  = diff(ip, x, y, z, 1, 0, 0) 
diff_y(ip, x, y, z)  = diff(ip, x, y, z, 0, 1, 0) 
diff_z(ip, x, y, z)  = diff(ip, x, y, z, 0, 0, 1) 

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

type dxdydz{M, N, O}
end

#=
@generated function poly{T<:Number, M, N, O}(ks::NTuple{64 , T}, 
                        x::T, y::T, z::T, 
                        trait::dxdydz{M, N, O}=dxdydz{0, 0, 0}())
    xt  = xterms[M+1]
    yt  = yterms[N+1]
    zt  = zterms[O+1]

    e   = :0
    sv  = "begin\n"
    for i = 0:3
        for j = 0:3
            sv  = "$sv x$(i)y$(j) = xt$i*yt$j\n" #mixted terms
            xiyj = parse("x$(i)y$(j)")
            for k = 0:3
                zk  = parse("zt$k")
                e   = :($xiyj*$zk*ks[$i*16 + $j*4 + $k + 1] + $e)
            end
        end
    end
    esv = parse("$sv end")
    return :($xt; $yt; $zt; $esv; $e)
end
=#

@generated function poly{T<:Number, M, N, O}(ks::NTuple{64, T}, 
                        x::T, y::T, z::T, 
                        trait::dxdydz{M, N, O}=dxdydz{0, 0, 0}())
    xt  = xterms[M+1]
    yt  = yterms[N+1]
    zt  = zterms[O+1]

    e   = :0
    for i = 0:3
        xi  = parse("xt$i")
        for j = 0:3
            yj  = parse("yt$j")
            for k = 0:3
                zk  = parse("zt$k")
                e   = :($xi*$yj*$zk*ks[$i*16 + $j*4 + $k + 1] + $e)
            end
        end
    end
    return :($xt; $yt; $zt; $e)
end

#coefficients for one cell
function cellcoeffs{T<:Number}(fs::Array{T, 3}, 
                    fs_x::Array{T, 3}, fs_y::Array{T, 3}, fs_z::Array{T, 3}, 
                    fs_xy::Array{T, 3}, fs_xz::Array{T, 3}, fs_yz::Array{T, 3},
                    fs_xyz::Array{T, 3})
    #code generated by maple
    t1 = 2*fs_z[1, 1, 1]
    t2 = 3*fs[1, 1, 1]
    t5 = 2*fs[1, 1, 1]
    t8 = 2*fs_yz[1, 1, 1]
    t9 = 3*fs_y[1, 1, 1]
    t10 = 3*fs_y[1, 1, 2]
    t12 = 2*fs_y[1, 1, 1]
    t13 = 2*fs_y[1, 1, 2]
    t17 = 3*fs_z[1, 1, 1]
    t18 = 3*fs_z[1, 2, 1]
    t20 = 6*fs_y[1, 1, 1]
    t21 = 6*fs_y[1, 1, 2]
    t22 = 3*fs_y[1, 2, 1]
    t23 = 3*fs_y[1, 2, 2]
    t24 = 4*fs_yz[1, 1, 1]
    t25 = 2*fs_yz[1, 1, 2]
    t26 = 2*fs_yz[1, 2, 1]
    t27 = 6*fs_z[1, 1, 1]
    t28 = 3*fs_z[1, 1, 2]
    t29 = 6*fs_z[1, 2, 1]
    t30 = 3*fs_z[1, 2, 2]
    t31 = 9*fs[1, 1, 1]
    t32 = 9*fs[1, 1, 2]
    t33 = 9*fs[1, 2, 1]
    t35 = t20 - t21 + t22 - t23 + t24 + t25 + t26 + fs_yz[1, 2, 2] + t27 + t28 - t29 - t30 + t31 - t32 - t33 + 9*fs[1, 2, 2]
    t36 = 4*fs_y[1, 1, 1]
    t37 = 4*fs_y[1, 1, 2]
    t38 = 6*fs[1, 1, 1]
    t39 = 6*fs[1, 1, 2]
    t40 = 6*fs[1, 2, 1]
    t41 = 6*fs[1, 2, 2]
    t42 = 2*fs_y[1, 2, 2]
    t43 = 2*fs_y[1, 2, 1]
    t44 = -t8 - t36 + t37 - t38 + t39 + t40 - t17 + t18 - t41 - fs_yz[1, 2, 1] + t42 - t43 + t30 - t25 - t28 - fs_yz[1, 2, 2]
    t47 = 2*fs_z[1, 2, 1]
    t49 = 4*fs_z[1, 1, 1]
    t50 = 4*fs_z[1, 2, 1]
    t51 = 2*fs_z[1, 2, 2]
    t52 = 2*fs_z[1, 1, 2]
    t53 = -t49 - t8 + t50 - t38 + t39 - t9 + t10 + t40 - t41 - t26 + t23 - t22 + t51 - fs_yz[1, 1, 2] - t52 - fs_yz[1, 2, 2]
    t54 = 4*fs[1, 1, 1]
    t55 = 4*fs[1, 1, 2]
    t56 = 4*fs[1, 2, 1]
    t58 = t54 - t55 + t12 + fs_yz[1, 1, 1] - t13 - t56 + t1 - t47 + 4*fs[1, 2, 2] + fs_yz[1, 2, 1] - t42 + t43 - t51 + fs_yz[1, 1, 2] + t52 + fs_yz[1, 2, 2]
    t59 = 2*fs_xz[1, 1, 1]
    t60 = 3*fs_x[1, 1, 1]
    t61 = 3*fs_x[1, 1, 2]
    t63 = 2*fs_x[1, 1, 1]
    t64 = 2*fs_x[1, 1, 2]
    t66 = 2*fs_xyz[1, 1, 1]
    t67 = 3*fs_xy[1, 1, 1]
    t68 = 3*fs_xy[1, 1, 2]
    t70 = 2*fs_xy[1, 1, 1]
    t71 = 2*fs_xy[1, 1, 2]
    t73 = 3*fs_x[1, 2, 1]
    t75 = 3*fs_xz[1, 1, 1]
    t76 = 3*fs_xz[1, 2, 1]
    t78 = 6*fs_xz[1, 1, 1]
    t79 = 3*fs_xz[1, 1, 2]
    t80 = 6*fs_xz[1, 2, 1]
    t81 = 3*fs_xz[1, 2, 2]
    t86 = 6*fs_xy[1, 1, 1]
    t87 = 6*fs_xy[1, 1, 2]
    t88 = 3*fs_xy[1, 2, 1]
    t89 = 3*fs_xy[1, 2, 2]
    t90 = 4*fs_xyz[1, 1, 1]
    t91 = 2*fs_xyz[1, 1, 2]
    t92 = 2*fs_xyz[1, 2, 1]
    t93 = t78 + t79 - t80 - t81 + 9*fs_x[1, 1, 1] - 9*fs_x[1, 1, 2] - 9*fs_x[1, 2, 1] + 9*fs_x[1, 2, 2] + t86 - t87 + t88 - t89 + t90 + t91 + t92 + fs_xyz[1, 2, 2]
    t94 = 4*fs_xy[1, 1, 1]
    t95 = 4*fs_xy[1, 1, 2]
    t96 = 6*fs_x[1, 1, 1]
    t97 = 6*fs_x[1, 1, 2]
    t98 = 6*fs_x[1, 2, 1]
    t99 = 6*fs_x[1, 2, 2]
    t100 = 2*fs_xy[1, 2, 2]
    t101 = 2*fs_xy[1, 2, 1]
    t102 = -t66 - t94 + t95 - t96 + t97 + t98 - t75 + t76 - t99 - fs_xyz[1, 2, 1] + t100 - t101 + t81 - t91 - t79 - fs_xyz[1, 2, 2]
    t103 = 2*fs_x[1, 2, 1]
    t105 = 2*fs_xz[1, 2, 1]
    t107 = 4*fs_xz[1, 1, 1]
    t108 = 4*fs_xz[1, 2, 1]
    t109 = 2*fs_xz[1, 2, 2]
    t110 = 2*fs_xz[1, 1, 2]
    t111 = -t107 - t66 + t108 - t96 + t97 - t67 + t68 + t98 - t99 - t92 + t89 - t88 + t109 - fs_xyz[1, 1, 2] - t110 - fs_xyz[1, 2, 2]
    t112 = 4*fs_x[1, 1, 1]
    t113 = 4*fs_x[1, 1, 2]
    t114 = 4*fs_x[1, 2, 1]
    t116 = t112 - t113 + t70 + fs_xyz[1, 1, 1] - t71 - t114 + t59 - t105 + 4*fs_x[1, 2, 2] + fs_xyz[1, 2, 1] - t100 + t101 - t109 + fs_xyz[1, 1, 2] + t110 + fs_xyz[1, 2, 2]
    t119 = 3*fs_z[2, 1, 1]
    t121 = 2*fs_xz[2, 1, 1]
    t122 = 6*fs_z[2, 1, 1]
    t123 = 3*fs_z[2, 1, 2]
    t124 = 9*fs[2, 1, 1]
    t126 = 3*fs_x[2, 1, 1]
    t127 = 3*fs_x[2, 1, 2]
    t128 = t107 + t110 + t121 + fs_xz[2, 1, 2] + t27 + t28 - t122 - t123 + t31 - t32 - t124 + 9*fs[2, 1, 2] + t96 - t97 + t126 - t127
    t129 = 6*fs[2, 1, 1]
    t130 = 6*fs[2, 1, 2]
    t131 = 2*fs_x[2, 1, 2]
    t132 = 2*fs_x[2, 1, 1]
    t133 = -t59 - t112 + t113 - t38 + t39 + t129 - t17 + t119 - t130 - fs_xz[2, 1, 1] + t131 - t132 + t123 - t110 - t28 - fs_xz[2, 1, 2]
    t134 = 3*fs_y[2, 1, 1]
    t136 = 3*fs_yz[1, 1, 1]
    t137 = 3*fs_yz[2, 1, 1]
    t139 = 2*fs_xyz[2, 1, 1]
    t140 = 9*fs_y[1, 1, 1]
    t141 = 9*fs_y[1, 1, 2]
    t142 = 9*fs_y[2, 1, 1]
    t143 = 9*fs_y[2, 1, 2]
    t144 = 6*fs_yz[1, 1, 1]
    t145 = 3*fs_yz[1, 1, 2]
    t146 = 6*fs_yz[2, 1, 1]
    t147 = 3*fs_yz[2, 1, 2]
    t148 = 3*fs_xy[2, 1, 1]
    t149 = 3*fs_xy[2, 1, 2]
    t150 = t139 + fs_xyz[2, 1, 2] + t140 - t141 - t142 + t143 + t144 + t145 - t146 - t147 + t86 - t87 + t148 - t149 + t90 + t91
    t151 = 6*fs_y[2, 1, 1]
    t152 = 6*fs_y[2, 1, 2]
    t153 = 2*fs_xy[2, 1, 2]
    t154 = 2*fs_xy[2, 1, 1]
    t155 = -t66 - t94 + t95 - t20 + t21 + t151 - t136 + t137 - t152 - fs_xyz[2, 1, 1] + t153 - t154 + t147 - t91 - t145 - fs_xyz[2, 1, 2]
    t156 = 3*fs_y[2, 2, 1]
    t158 = 3*fs_x[2, 2, 1]
    t159 = t20 + t22 - t151 - t156 + t31 - t33 - t124 + 9*fs[2, 2, 1] + t96 - t98 + t126 - t158 + t94 + t101 + t154 + fs_xy[2, 2, 1]
    t160 = 9*fs_z[2, 2, 1]
    t161 = 3*fs_xz[2, 1, 1]
    t162 = 3*fs_xz[2, 2, 1]
    t163 = 3*fs_yz[1, 2, 1]
    t164 = 3*fs_yz[2, 2, 1]
    t165 = 9*fs_z[1, 1, 1]
    t166 = 9*fs_z[1, 2, 1]
    t167 = 9*fs_z[2, 1, 1]
    t168 = t160 + t139 + fs_xyz[2, 2, 1] + t78 - t80 + t161 - t162 + t144 + t163 - t146 - t164 + t165 - t166 - t167 + t90 + t92
    t169 = 6*fs_xz[2, 2, 1]
    t170 = 6*fs_yz[2, 2, 1]
    t172 = 9*fs_x[2, 2, 2]
    t173 = 9*fs_y[2, 2, 2]
    t174 = 9*fs_z[2, 2, 2]
    t175 = 3*fs_xy[2, 2, 2]
    t176 = 3*fs_xz[2, 2, 2]
    t177 = 3*fs_yz[2, 2, 2]
    t178 = 9*fs_y[2, 2, 1]
    t180 = 9*fs_x[2, 2, 1]
    t181 = 3*fs_xy[2, 2, 1]
    t183 = 2*fs_xyz[2, 2, 1]
    t184 = t169 + t170 + 27*fs[2, 2, 2] - t172 - t173 - t174 + t175 + t176 + t177 - fs_xyz[2, 2, 2] + t178 - 27*fs[2, 2, 1] + t180 - t181 - 18*fs_z[2, 2, 1] - t183
    t186 = 6*fs_xy[2, 1, 1]
    t188 = 4*fs_xyz[2, 1, 1]
    t189 = 2*fs_xyz[2, 1, 2]
    t191 = 6*fs_yz[2, 1, 2]
    t192 = 6*fs_xy[2, 1, 2]
    t194 = 6*fs_xz[2, 1, 1]
    t195 = 3*fs_xz[2, 1, 2]
    t196 = 9*fs_z[2, 1, 2]
    t198 = 9*fs_x[2, 1, 2]
    t200 = 9*fs_x[2, 1, 1]
    t201 = 18*fs_y[2, 1, 1] - t186 + 12*fs_yz[2, 1, 1] - t188 - t189 - 18*fs_y[2, 1, 2] + t191 + t192 + 18*fs_z[2, 1, 1] - t194 - t195 + t196 - 27*fs[2, 1, 2] + t198 + 27*fs[2, 1, 1] - t200
    t203 = 6*fs_xy[1, 2, 2]
    t204 = 2*fs_xyz[1, 2, 2]
    t206 = 4*fs_xyz[1, 1, 2]
    t208 = 6*fs_xy[1, 2, 1]
    t210 = 4*fs_xyz[1, 2, 1]
    t211 = 6*fs_xz[1, 2, 2]
    t216 = 6*fs_xz[1, 1, 2]
    t219 = t203 - t204 + 12*fs_xy[1, 1, 2] - t206 + 18*fs_x[1, 2, 1] - t208 + 12*fs_xz[1, 2, 1] - t210 + t211 - 18*fs_x[1, 2, 2] - 18*fs_x[1, 1, 1] - 12*fs_xz[1, 1, 1] + 18*fs_x[1, 1, 2] - t216 - 12*fs_xy[1, 1, 1] - 8*fs_xyz[1, 1, 1]
    t221 = 9*fs_y[1, 2, 1]
    t223 = 6*fs_yz[1, 2, 1]
    t224 = 9*fs_y[1, 2, 2]
    t225 = 3*fs_yz[1, 2, 2]
    t226 = 9*fs_z[1, 2, 2]
    t231 = 9*fs_z[1, 1, 2]
    t235 = 6*fs_yz[1, 1, 2]
    t236 = 27*fs[1, 2, 1] - t221 + 18*fs_z[1, 2, 1] - t223 + t224 - t225 + t226 - 27*fs[1, 2, 2] - 27*fs[1, 1, 1] - 18*fs_z[1, 1, 1] + 27*fs[1, 1, 2] - t231 - 18*fs_y[1, 1, 1] - 12*fs_yz[1, 1, 1] + 18*fs_y[1, 1, 2] - t235
    t239 = 18*fs[2, 2, 2]
    t240 = 6*fs_x[2, 2, 2]
    t241 = 6*fs_y[2, 2, 2]
    t242 = 2*fs_xy[2, 2, 2]
    t243 = 6*fs_y[2, 2, 1]
    t244 = 18*fs[2, 2, 1]
    t245 = 6*fs_x[2, 2, 1]
    t246 = 2*fs_xy[2, 2, 1]
    t247 = -t162 - t164 - t239 + t240 + t241 + t174 - t242 - t176 - t177 + fs_xyz[2, 2, 2] - t243 + t244 - t245 + t246 + t160 + fs_xyz[2, 2, 1]
    t248 = 12*fs_y[2, 1, 1]
    t249 = 4*fs_xy[2, 1, 1]
    t250 = 12*fs_y[2, 1, 2]
    t251 = 4*fs_xy[2, 1, 2]
    t252 = 18*fs[2, 1, 2]
    t253 = 6*fs_x[2, 1, 2]
    t254 = 18*fs[2, 1, 1]
    t255 = 6*fs_x[2, 1, 1]
    t256 = -t248 + t249 - t146 + t139 + t189 + t250 - t191 - t251 - t167 + t161 + t195 - t196 + t252 - t253 - t254 + t255
    t258 = 4*fs_xy[1, 2, 2]
    t260 = 12*fs_x[1, 2, 1]
    t261 = 4*fs_xy[1, 2, 1]
    t262 = 12*fs_x[1, 2, 2]
    t263 = 12*fs_x[1, 1, 1]
    t264 = 12*fs_x[1, 1, 2]
    t266 = -t258 + t204 - 8*fs_xy[1, 1, 2] + t206 - t260 + t261 - t80 + t92 - t211 + t262 + t263 + t78 - t264 + t216 + 8*fs_xy[1, 1, 1] + t90
    t267 = 18*fs[1, 2, 1]
    t268 = 6*fs_y[1, 2, 1]
    t269 = 6*fs_y[1, 2, 2]
    t270 = 18*fs[1, 2, 2]
    t271 = 18*fs[1, 1, 1]
    t272 = 18*fs[1, 1, 2]
    t273 = 12*fs_y[1, 1, 1]
    t274 = 12*fs_y[1, 1, 2]
    t275 = -t267 + t268 - t166 + t163 - t269 + t225 - t226 + t270 + t271 + t165 - t272 + t231 + t273 + t144 - t274 + t235
    t278 = 6*fs[2, 2, 1]
    t279 = 2*fs_x[2, 2, 1]
    t280 = t129 - t9 - t70 + t134 - t38 + t40 - t112 + t114 - t278 - fs_xy[2, 1, 1] + t279 - t132 + t156 - t101 - t22 - fs_xy[2, 2, 1]
    t281 = 6*fs_z[2, 2, 1]
    t282 = 2*fs_xz[2, 2, 1]
    t283 = -t107 + t122 - t136 - t66 + t137 - t27 + t29 + t108 - t281 - fs_xyz[2, 1, 1] + t282 - t121 + t164 - t92 - t163 - fs_xyz[2, 2, 1]
    t284 = 4*fs_xz[2, 2, 1]
    t285 = 6*fs_z[2, 2, 2]
    t286 = 2*fs_xz[2, 2, 2]
    t287 = 12*fs_z[2, 2, 1]
    t288 = -t284 - t170 - t239 + t240 + t173 + t285 - t175 - t286 - t177 + fs_xyz[2, 2, 2] - t178 + t244 - t245 + t181 + t287 + t183
    t289 = 12*fs_z[2, 1, 1]
    t290 = 4*fs_xz[2, 1, 1]
    t291 = 2*fs_xz[2, 1, 2]
    t292 = 6*fs_z[2, 1, 2]
    t293 = -t142 + t148 - t146 + t139 + fs_xyz[2, 1, 2] + t143 - t147 - t149 - t289 + t290 + t291 - t292 + t252 - t253 - t254 + t255
    t296 = 4*fs_xz[1, 2, 2]
    t298 = 4*fs_xz[1, 1, 2]
    t299 = -t203 + t204 - t87 + t91 - t260 + t208 - 8*fs_xz[1, 2, 1] + t210 - t296 + t262 + t263 + 8*fs_xz[1, 1, 1] - t264 + t298 + t86 + t90
    t300 = 12*fs_z[1, 2, 1]
    t301 = 6*fs_z[1, 2, 2]
    t302 = 12*fs_z[1, 1, 1]
    t303 = 6*fs_z[1, 1, 2]
    t304 = -t267 + t221 - t300 + t223 - t224 + t225 - t301 + t270 + t271 + t302 - t272 + t303 + t140 + t144 - t141 + t145
    t307 = 12*fs[2, 2, 2]
    t308 = 4*fs_x[2, 2, 2]
    t309 = 12*fs[2, 2, 1]
    t310 = 4*fs_x[2, 2, 1]
    t311 = t282 + t164 + t307 - t308 - t241 - t285 + t242 + t286 + t177 - fs_xyz[2, 2, 2] + t243 - t309 + t310 - t246 - t281 - fs_xyz[2, 2, 1]
    t312 = 12*fs[2, 1, 2]
    t313 = 4*fs_x[2, 1, 2]
    t314 = 12*fs[2, 1, 1]
    t315 = 4*fs_x[2, 1, 1]
    t316 = t151 - t154 + t137 - fs_xyz[2, 1, 1] - fs_xyz[2, 1, 2] - t152 + t147 + t153 + t122 - t121 - t291 + t292 - t312 + t313 + t314 - t315
    t322 = t258 - t204 + t95 - t91 + 8*fs_x[1, 2, 1] - t261 + t108 - t92 + t296 - 8*fs_x[1, 2, 2] - 8*fs_x[1, 1, 1] - t107 + 8*fs_x[1, 1, 2] - t298 - t94 - t66
    t323 = 12*fs[1, 2, 1]
    t324 = 12*fs[1, 2, 2]
    t325 = 12*fs[1, 1, 1]
    t326 = 12*fs[1, 1, 2]
    t327 = t323 - t268 + t29 - t163 + t269 - t225 + t301 - t324 - t325 - t27 + t326 - t303 - t20 - t136 + t21 - t145
    t332 = 2*fs_z[2, 1, 1]
    t334 = 4*fs_z[2, 1, 1]
    t335 = 2*fs_z[2, 1, 2]
    t336 = -t49 - t59 + t334 - t38 + t39 - t60 + t61 + t129 - t130 - t121 + t127 - t126 + t335 - fs_xz[1, 1, 2] - t52 - fs_xz[2, 1, 2]
    t337 = 4*fs[2, 1, 1]
    t339 = t54 - t55 + t63 + fs_xz[1, 1, 1] - t64 - t337 + t1 - t332 + 4*fs[2, 1, 2] + fs_xz[2, 1, 1] - t131 + t132 - t335 + fs_xz[1, 1, 2] + t52 + fs_xz[2, 1, 2]
    t340 = 2*fs_y[2, 1, 1]
    t342 = 2*fs_yz[2, 1, 1]
    t344 = 4*fs_yz[2, 1, 1]
    t345 = 2*fs_yz[2, 1, 2]
    t346 = -t24 - t66 + t344 - t20 + t21 - t67 + t68 + t151 - t152 - t139 + t149 - t148 + t345 - fs_xyz[1, 1, 2] - t25 - fs_xyz[2, 1, 2]
    t347 = 4*fs_y[2, 1, 1]
    t348 = 4*fs_y[2, 1, 2]
    t349 = t36 - t37 + t70 + fs_xyz[1, 1, 1] - t71 - t347 + t8 - t342 + t348 + fs_xyz[2, 1, 1] - t153 + t154 - t345 + fs_xyz[1, 1, 2] + t25 + fs_xyz[2, 1, 2]
    t350 = 2*fs_y[2, 2, 1]
    t351 = -t36 - t70 + t347 + t129 - t38 + t40 - t60 + t73 - t278 - t154 + t158 - t126 + t350 - fs_xy[1, 2, 1] - t43 - fs_xy[2, 2, 1]
    t352 = 2*fs_yz[2, 2, 1]
    t353 = -t24 - t66 + t344 - t75 + t122 - t27 + t29 + t76 - t281 - t139 + t162 - t161 + t352 - fs_xyz[1, 2, 1] - t26 - fs_xyz[2, 2, 1]
    t354 = 4*fs_yz[2, 2, 1]
    t355 = 2*fs_yz[2, 2, 2]
    t356 = -t169 - t354 - t239 + t172 + t241 + t285 - t175 - t176 - t355 + fs_xyz[2, 2, 2] - t243 + t244 - t180 + t181 + t287 + t183
    t358 = 4*fs_yz[2, 1, 2]
    t359 = -t248 + t186 - 8*fs_yz[2, 1, 1] + t188 + t189 + t250 - t358 - t192 - t289 + t194 + t195 - t292 + t252 - t198 - t254 + t200
    t361 = 4*fs_yz[1, 2, 1]
    t362 = 2*fs_yz[1, 2, 2]
    t364 = 4*fs_yz[1, 1, 2]
    t365 = -t267 + t268 - t300 + t361 - t269 + t362 - t301 + t270 + t271 + t302 - t272 + t303 + t273 + 8*fs_yz[1, 1, 1] - t274 + t364
    t368 = 4*fs_y[2, 2, 2]
    t369 = 4*fs_y[2, 2, 1]
    t370 = t162 + t352 + t307 - t240 - t368 - t285 + t242 + t176 + t355 - fs_xyz[2, 2, 2] + t369 - t309 + t245 - t246 - t281 - fs_xyz[2, 2, 1]
    t373 = 8*fs_y[2, 1, 1] - t249 + t344 - t139 - t189 - 8*fs_y[2, 1, 2] + t358 + t251 + t122 - t161 - t195 + t292 - t312 + t253 + t314 - t255
    t375 = 4*fs_y[1, 2, 1]
    t376 = 4*fs_y[1, 2, 2]
    t379 = t323 - t375 + t29 - t26 + t376 - t362 + t301 - t324 - t325 - t27 + t326 - t303 - 8*fs_y[1, 1, 1] - t24 + 8*fs_y[1, 1, 2] - t364
    t383 = -t337 + t12 + fs_xy[1, 1, 1] - t340 + t54 - t56 + t63 - t103 + 4*fs[2, 2, 1] + fs_xy[2, 1, 1] - t279 + t132 - t350 + fs_xy[1, 2, 1] + t43 + fs_xy[2, 2, 1]
    t384 = 4*fs_z[2, 2, 1]
    t385 = t59 - t334 + t8 + fs_xyz[1, 1, 1] - t342 + t49 - t50 - t105 + t384 + fs_xyz[2, 1, 1] - t282 + t121 - t352 + fs_xyz[1, 2, 1] + t26 + fs_xyz[2, 2, 1]
    t386 = 4*fs_z[2, 2, 2]
    t388 = t284 + t354 + t307 - t240 - t241 - t386 + t175 + t286 + t355 - fs_xyz[2, 2, 2] + t243 - t309 + t245 - t181 - 8*fs_z[2, 2, 1] - t183
    t390 = 4*fs_z[2, 1, 2]
    t391 = t151 - t148 + t344 - t139 - fs_xyz[2, 1, 2] - t152 + t345 + t149 + 8*fs_z[2, 1, 1] - t290 - t291 + t390 - t312 + t253 + t314 - t255
    t394 = 4*fs_z[1, 2, 2]
    t396 = 4*fs_z[1, 1, 2]
    t397 = t323 - t268 + 8*fs_z[1, 2, 1] - t361 + t269 - t362 + t394 - t324 - t325 - 8*fs_z[1, 1, 1] + t326 - t396 - t20 - t24 + t21 - t25
    t402 = -t282 - t352 - 8*fs[2, 2, 2] + t308 + t368 + t386 - t242 - t286 - t355 + fs_xyz[2, 2, 2] - t369 + 8*fs[2, 2, 1] - t310 + t246 + t384 + fs_xyz[2, 2, 1]
    t405 = -t347 + t154 - t342 + fs_xyz[2, 1, 1] + fs_xyz[2, 1, 2] + t348 - t345 - t153 - t334 + t121 + t291 - t390 + 8*fs[2, 1, 2] - t313 - 8*fs[2, 1, 1] + t315
    t411 = -8*fs[1, 2, 1] + t375 - t50 + t26 - t376 + t362 - t394 + 8*fs[1, 2, 2] + 8*fs[1, 1, 1] + t49 - 8*fs[1, 1, 2] + t396 + t36 + t8 - t37 + t25

    return (
        fs[1, 1, 1],
        fs_z[1, 1, 1],
        -t1 - t2 + 3*fs[1, 1, 2] - fs_z[1, 1, 2],
        t5 + fs_z[1, 1, 1] - 2*fs[1, 1, 2] + fs_z[1, 1, 2],
        fs_y[1, 1, 1],
        fs_yz[1, 1, 1],
        -t8 - t9 + t10 - fs_yz[1, 1, 2],
        t12 + fs_yz[1, 1, 1] - t13 + fs_yz[1, 1, 2],
        -t12 - t2 + 3*fs[1, 2, 1] - fs_y[1, 2, 1],
        -t8 - t17 + t18 - fs_yz[1, 2, 1],
        t35,
        t44,
        t5 + fs_y[1, 1, 1] - 2*fs[1, 2, 1] + fs_y[1, 2, 1],
        t1 + fs_yz[1, 1, 1] - t47 + fs_yz[1, 2, 1],
        t53,
        t58,
        fs_x[1, 1, 1],
        fs_xz[1, 1, 1],
        -t59 - t60 + t61 - fs_xz[1, 1, 2],
        t63 + fs_xz[1, 1, 1] - t64 + fs_xz[1, 1, 2],
        fs_xy[1, 1, 1],
        fs_xyz[1, 1, 1],
        -t66 - t67 + t68 - fs_xyz[1, 1, 2],
        t70 + fs_xyz[1, 1, 1] - t71 + fs_xyz[1, 1, 2],
        -t70 - t60 + t73 - fs_xy[1, 2, 1],
        -t66 - t75 + t76 - fs_xyz[1, 2, 1],
        t93,
        t102,
        t63 + fs_xy[1, 1, 1] - t103 + fs_xy[1, 2, 1],
        t59 + fs_xyz[1, 1, 1] - t105 + fs_xyz[1, 2, 1],
        t111,
        t116,
        -t63 - t2 + 3*fs[2, 1, 1] - fs_x[2, 1, 1],
        -t59 - t17 + t119 - fs_xz[2, 1, 1],
        t128,
        t133,
        -t70 - t9 + t134 - fs_xy[2, 1, 1],
        -t66 - t136 + t137 - fs_xyz[2, 1, 1],
        t150,
        t155,
        t159,
        t168,
        t184 + t201 + t219 + t236,
        t247 + t256 + t266 + t275,
        t280,
        t283,
        t288 + t293 + t299 + t304,
        t311 + t316 + t322 + t327,
        t5 + fs_x[1, 1, 1] - 2*fs[2, 1, 1] + fs_x[2, 1, 1],
        t1 + fs_xz[1, 1, 1] - t332 + fs_xz[2, 1, 1],
        t336,
        t339,
        t12 + fs_xy[1, 1, 1] - t340 + fs_xy[2, 1, 1],
        t8 + fs_xyz[1, 1, 1] - t342 + fs_xyz[2, 1, 1],
        t346,
        t349,
        t351,
        t353,
        t356 + t359 + t93 + t365,
        t370 + t373 + t102 + t379,
        t383,
        t385,
        t388 + t391 + t111 + t397,
        t402 + t405 + t116 + t411,
    )
end

