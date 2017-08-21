################################################################################
#  Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

function basic4d()
    tols= [1e-3, 1e-2, 1e-1, 1]

    #Test functions
    A   = 2
    B   = 1//2
    C   = 3//4
    D   = 2//5

    funs = [ 
        (x, y, z, w) -> x + y + z + w,
        (x, y, z, w) -> x*y*z*w,
        (x, y, z, w) -> x^2*y^2*z^2*w^2,
        (x, y, z, w) -> log(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1) 
    ]

    funs_x = [
        (x, y, z, w) -> one(x),
        (x, y, z, w) -> y*z*w,
        (x, y, z, w) -> 2x*y^2*z^2*w^2,
        (x, y, z, w) -> 2*A*x/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)
    ]

    funs_y = [
        (x, y, z, w) -> one(x),
        (x, y, z, w) -> x*z*w,
        (x, y, z, w) -> x^2*2y*z^2*w^2,
        (x, y, z, w) -> 2*B*y/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)
    ]

    funs_z = [
        (x, y, z, w) -> one(x),
        (x, y, z, w) -> x*y*w,
        (x, y, z, w) -> x^2*y^2*2z*w^2,
        (x, y, z, w) -> 2*C*z/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)
    ]

    funs_w = [
        (x, y, z, w) -> one(x),
        (x, y, z, w) -> x*y*z,
        (x, y, z, w) -> x^2*y^2*z^2*2w,
        (x, y, z, w) -> 2*D*w/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)
    ]

    funs_xy = [
        (x, y, z, w) -> zero(x),
        (x, y, z, w) -> z*w,
        (x, y, z, w) -> 4*x*y*z^2*w^2,
        (x, y, z, w) -> -4*A*B*x*y/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)^2
    ]

    funs_xz = [
        (x, y, z, w) -> zero(x),
        (x, y, z, w) -> y*w,
        (x, y, z, w) -> 4*x*y^2*z*w^2,
        (x, y, z, w) -> -4*A*C*x*z/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)^2
    ]

    funs_xw = [
        (x, y, z, w) -> zero(x),
        (x, y, z, w) -> y*z,
        (x, y, z, w) -> 4*x*y^2*z^2*w,
        (x, y, z, w) -> -4*A*D*x*w/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)^2
    ]

    funs_yz = [
        (x, y, z, w) -> zero(x),
        (x, y, z, w) -> x*w,
        (x, y, z, w) -> 4*x^2*y*z*w^2,
        (x, y, z, w) -> -4*B*C*y*z/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)^2
    ]

    funs_yw = [
        (x, y, z, w) -> zero(x),
        (x, y, z, w) -> x*z,
        (x, y, z, w) -> 4*x^2*y*z^2*w,
        (x, y, z, w) -> -4*B*D*y*w/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)^2
    ]

    funs_zw = [
        (x, y, z, w) -> zero(x),
        (x, y, z, w) -> x*y,
        (x, y, z, w) -> 4*x^2*y^2*z*w,
        (x, y, z, w) -> -4*C*D*z*w/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)^2
    ]

    funs_xyz = [
        (x, y, z, w) -> zero(x),
        (x, y, z, w) -> w,
        (x, y, z, w) -> 8*x*y*z*w^2,
        (x, y, z, w) -> 16*A*B*C*x*y*z/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)^3
    ]

    funs_xyw = [
        (x, y, z, w) -> zero(x),
        (x, y, z, w) -> z,
        (x, y, z, w) -> 8*x*y*z^2*w,
        (x, y, z, w) -> 16*A*B*D*x*y*w/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)^3
    ]

    funs_xzw = [
        (x, y, z, w) -> zero(x),
        (x, y, z, w) -> y,
        (x, y, z, w) -> 8*x*y^2*z*w,
        (x, y, z, w) -> 16*A*C*D*x*z*w/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)^3
    ]

    funs_yzw = [
        (x, y, z, w) -> zero(x),
        (x, y, z, w) -> x,
        (x, y, z, w) -> 8*x^2*y*z*w,
        (x, y, z, w) -> 16*B*C*D*y*z*w/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)^3
    ]

    funs_xyzw = [
        (x, y, z, w) -> zero(x),
        (x, y, z, w) -> one(x),
        (x, y, z, w) -> 16*x*y*z*w,
        (x, y, z, w) -> 
               -96*A*B*C*D*x*y*z*w/(A*x^2 + B*y^2 + C*z^2 + D*w^2 + 1)^4
    ]


    @testset "4D Basic" begin
        @out "\n4D Basic"
        @testset "$T" for T in [Float32 Float64]
            @out "$T" 1

            #Test interval
            xs  = 1:T(0.5):2
            xs0 = xs[1:end-1]
            xst = xs0 + step(xs)/3

            ys  = 5:T(1.):7
            ys0 = ys[1:end-1]'
            yst = ys0 + step(ys)/3

            zs  = 2:T(2.):6
            zs0 = Array{T}(1, 1, length(zs) - 1)
            zs0[1, 1, :] = zs[1:end-1]
            zst = zs0 + step(zs)/3

            ws  = 3:T(1.5):6
            ws0 = Array{T}(1, 1, 1, length(ws) - 1)
            ws0[1, 1, 1, :] = ws[1:end-1]
            wst = ws0 + step(ws)/3
            
            @testset "$b" for b in instances(Boundary)
                @out "$b" 2


                @testset "$f" for (f, f_x, f_y, f_z, f_w, f_xy, f_xz, f_xw, 
                                   f_yz, f_yw, f_zw, f_xyz, f_xyw, f_xzw,
                                   f_yzw, f_xyzw) in 
                            zip(funs, funs_x, funs_y, funs_z, funs_w,
                                funs_xy, funs_xz, funs_xw, funs_yz, funs_yw,
                                funs_zw, funs_xyz, funs_xyw, funs_xzw, 
                                funs_yzw, funs_xyzw)
                    @out "$f" 3

                    ip          = Interpolation(xs, ys, zs, ws, f,
                                    f_x, f_y, f_z, f_w,
                                    f_xy, f_xz, f_xw, f_yz, f_yw, f_zw,
                                    f_xyz, f_xyw, f_xzw, f_yzw, f_xyzw, b)

                    ip(xs0[1], ys0[1], zs0[1], ws0[1])
                    #ipf(x, y, z, w)    = interpolate(ip, x, y, z, w)

                   ip_x(x, y, z, w)   = diff_x(ip, x, y, z, w)
                   ip_y(x, y, z, w)   = diff_y(ip, x, y, z, w)
                   ip_z(x, y, z, w)   = diff_z(ip, x, y, z, w)
                   ip_w(x, y, z, w)   = diff_w(ip, x, y, z, w)

                   ip_xy(x, y, z, w)  = diff(ip, x, y, z, w, 1, 1, 0, 0)
                   ip_xz(x, y, z, w)  = diff(ip, x, y, z, w, 1, 0, 1, 0)
                   ip_xw(x, y, z, w)  = diff(ip, x, y, z, w, 1, 0, 0, 1)
                   ip_yz(x, y, z, w)  = diff(ip, x, y, z, w, 0, 1, 1, 0)
                   ip_yw(x, y, z, w)  = diff(ip, x, y, z, w, 0, 1, 0, 1)
                   ip_zw(x, y, z, w)  = diff(ip, x, y, z, w, 0, 0, 1, 1)

                   ip_xyz(x, y, z, w) = diff(ip, x, y, z, w, 1, 1, 1, 0)
                   ip_xyw(x, y, z, w) = diff(ip, x, y, z, w, 1, 1, 0, 1)
                   ip_xzw(x, y, z, w) = diff(ip, x, y, z, w, 1, 0, 1, 1)
                   ip_yzw(x, y, z, w) = diff(ip, x, y, z, w, 0, 1, 1, 1)

                   ip_xyzw(x, y, z, w)= diff(ip, x, y, z, w, 1, 1, 1, 1)

                    @testset "interpolation at xp = 0" begin
                        @test isapprox(ip.(xs0, ys0, zs0, ws0), 
                                f.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_x.(xs0, ys0, zs0, ws0), 
                                f_x.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_y.(xs0, ys0, zs0, ws0), 
                                f_y.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_z.(xs0, ys0, zs0, ws0), 
                                f_z.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_w.(xs0, ys0, zs0, ws0), 
                                f_w.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_xy.(xs0, ys0, zs0, ws0), 
                                f_xy.(xs0, ys0, zs0, ws0), atol=tols[3]) 
                        @test isapprox(ip_xz.(xs0, ys0, zs0, ws0), 
                                f_xz.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_xw.(xs0, ys0, zs0, ws0), 
                                f_xw.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_yz.(xs0, ys0, zs0, ws0), 
                                f_yz.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_yw.(xs0, ys0, zs0, ws0), 
                                f_yw.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_zw.(xs0, ys0, zs0, ws0), 
                                f_zw.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_xyz.(xs0, ys0, zs0, ws0), 
                                f_xyz.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_xyw.(xs0, ys0, zs0, ws0), 
                                f_xyw.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_xzw.(xs0, ys0, zs0, ws0), 
                                f_xzw.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_yzw.(xs0, ys0, zs0, ws0), 
                                f_yzw.(xs0, ys0, zs0, ws0), atol=tols[3])
                        @test isapprox(ip_xyzw.(xs0, ys0, zs0, ws0), 
                                f_xyzw.(xs0, ys0, zs0, ws0), atol=tols[3])
                    end

                    @testset "interpolation in cells" begin
                        @test isapprox(ip.(xst, yst, zst, wst), 
                                        f.(xst, yst, zst, wst), 
                                        atol=tols[1], rtol=tols[1])
                        @test isapprox(ip_x.(xst, yst, zst, wst), 
                                        f_x.(xst, yst, zst, wst), 
                                        atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_y.(xst, yst, zst, wst), 
                                       f_y.(xst, yst, zst, wst), 
                                       atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_z.(xst, yst, zst, wst), 
                                       f_z.(xst, yst, zst, wst), 
                                       atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_w.(xst, yst, zst, wst), 
                                       f_w.(xst, yst, zst, wst), 
                                       atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_xy.(xst, yst, zst, wst), 
                                       f_xy.(xst, yst, zst, wst), 
                                       atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_xz.(xst, yst, zst, wst), 
                                       f_xz.(xst, yst, zst, wst), 
                                       atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_xw.(xst, yst, zst, wst), 
                                       f_xw.(xst, yst, zst, wst), 
                                       atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_yz.(xst, yst, zst, wst), 
                                       f_yz.(xst, yst, zst, wst), 
                                        atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_yw.(xst, yst, zst, wst), 
                                       f_yw.(xst, yst, zst, wst),   
                                       atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_zw.(xst, yst, zst, wst), 
                                       f_zw.(xst, yst, zst, wst), 
                                       atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_xyz.(xst, yst, zst, wst), 
                                       f_xyz.(xst, yst, zst, wst), 
                                       atol=tols[3], rtol=tols[3])
                        @test isapprox(ip_xyw.(xst, yst, zst, wst), 
                                       f_xyw.(xst, yst, zst, wst), 
                                       atol=tols[3], rtol=tols[3])
                        @test isapprox(ip_xzw.(xst, yst, zst, wst), 
                                       f_xzw.(xst, yst, zst, wst), 
                                       atol=tols[3], rtol=tols[3])
                        @test isapprox(ip_yzw.(xst, yst, zst, wst), 
                                       f_yzw.(xst, yst, zst, wst), 
                                       atol=tols[3], rtol=tols[3])
                        @test isapprox(ip_xyzw.(xst, yst, zst, wst), 
                                       f_xyzw.(xst, yst, zst, wst), 
                                       atol=tols[3], rtol=tols[3])
                    end
                end
            end
        end
    end
end

function exp4d()
    tols= [1e-2, 1e-1, 1e0]

    #Test interval
    xs  = 1:0.5:2
    xs0 = xs[1:end-1]
    xst = xs0 + step(xs)/3

    ys  = 5:1.:7
    ys0 = ys[1:end-1]'
    yst = ys0 + step(ys)/3

    zs  = 2:0.75:3.5
    zs0 = Array{Float64}(1, 1, length(zs) - 1)
    zs0[1, 1, :] = zs[1:end-1]
    zst = zs0 + step(zs)/3

    ws  = 3:0.4:4
    ws0 = Array{Float64}(1, 1, 1, length(ws) - 1)
    ws0[1, 1, 1, :] = ws[1:end-1]
    wst = ws0 + step(ws)/3

    exp4d(x, y, z, w) = exp(x + y + z + w)

    ip  = Interpolation(xs, ys, zs, ws, ntuple(i->exp4d, 16)...)

    @testset "4D exp" begin
        @out "\n4D exp"
        @test isapprox(ip.(xs0, ys0, zs0, ws0), exp4d.(xs0, ys0, zs0, ws0))
        @test isapprox(ip.(xst, yst, zst, wst), exp4d.(xst, yst, zst, wst), 
                       rtol=tols[1], atol=tols[1])

        for i = 0:2
         for j = 0:2
          for k = 0:2
           for l = 0:2
                n   = i+j+k+l
                if n > 0
                    m   = maximum([i, j, k, l]) + 1
                    @out "d^$(n)f/dx^$i/dy^$j/dz^$k/dw^$l" 1
                    dip(x, y, z, w)  = diff(ip, x, y, z, w, i, j, k, l)
                    @test isapprox(dip.(xst, yst, zst, wst), 
                                   exp4d.(xst, yst, zst, wst), 
                                   rtol=tols[m], atol=tols[m])
                end
           end
          end
         end
        end
    end
end

function boundary4d()
    tols= [1e-3, 1e-2, 1e-1, 1]

    fr  = 89//97
    big = 1e6

    xs  = 1:0.2:2
    ys  = 5:0.4:7
    zs  = 2:1.:7
    ws  = 3:0.6:6

    dx  = xs[end] - xs[1]
    dy  = ys[end] - ys[1]
    dz  = zs[end] - zs[1]
    dw  = ws[end] - ws[1]

    exp4d(x, y, z, w) = exp(x + y + z + w)

    xts = [
        xs[end]  + 1/big,
        xs[1]    - 1/big,
        xs[end] +  dx*fr,
        xs[1] -  dx*fr,
        xs[end] + big,
        xs[1] - big
    ]

    yts = [
        ys[end]  + 1/big,
        ys[1]    - 1/big,
        ys[end] +  dy*fr,
        ys[1] -  dy*fr,
        ys[end] + big,
        ys[1] - big
    ]

    zts = [
        zs[end]  + 1/big,
        zs[1]    - 1/big,
        zs[end] +  dz*fr,
        zs[1] -  dz*fr,
        zs[end] + big,
        zs[1] - big
    ]

    wts = [
        ws[end]  + 1/big,
        ws[1]    - 1/big,
        ws[end] +  dw*fr,
        ws[1] -  dw*fr,
        ws[end] + big,
        ws[1] - big
    ]

    @testset "4D Boundary" begin
        @out  "\n4D Boundary"
        @testset "Mode throw_error" begin
            @out  "throw_error" 1
            ip  = Interpolation(xs, ys, zs, ws, ntuple(i->exp4d, 16)...,
                  throw_error)

            ip_x(x, y, z, w)  = diff_x(ip, x, y, z, w)
            ip_y(x, y, z, w)  = diff_y(ip, x, y, z, w)
            ip_z(x, y, z, w)  = diff_z(ip, x, y, z, w)
            ip_w(x, y, z, w)  = diff_w(ip, x, y, z, w)

            for (xt, yt, zt, wt) in 
              collect(Base.product(xts, yts, zts, wts))[2:end]  
                @test_throws BoundsError ip(xt, yt, zt, wt)
                @test_throws BoundsError ip_x(xt, yt, zt, wt)
                @test_throws BoundsError ip_y(xt, yt, zt, wt)
                @test_throws BoundsError ip_z(xt, yt, zt, wt)
                @test_throws BoundsError ip_w(xt, yt, zt, wt)
            end
        end

        #NO TESTING FOR UNSAFE!

        @testset "Mode undef" begin
            @out  "undef" 1
            ip  = Interpolation(xs, ys, zs, ws, ntuple(i->exp4d, 16)...,
                  undef)

            ip_x(x, y, z, w)  = diff_x(ip, x, y, z, w)
            ip_y(x, y, z, w)  = diff_y(ip, x, y, z, w)
            ip_z(x, y, z, w)  = diff_z(ip, x, y, z, w)
            ip_w(x, y, z, w)  = diff_w(ip, x, y, z, w)

            for (xt, yt, zt, wt) in 
              collect(Base.product(xts, yts, zts, wts))[2:end]  
                @test isnan(ip(xt, yt, zt, wt))
                @test isnan(ip_x(xt, yt, zt, wt))
                @test isnan(ip_y(xt, yt, zt, wt))
                @test isnan(ip_z(xt, yt, zt, wt))
                @test isnan(ip_w(xt, yt, zt, wt))
            end
        end
    end
end
