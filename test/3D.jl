################################################################################
#  Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

#Test sets
function basic3d()
    tols= [1e-3, 1e-2, 1e-1, 1]

    #Test functions
    A   = 2
    B   = 1//2
    C   = 3//4

    funs = [ 
        (x, y, z) -> x + y + z,
        (x, y, z) -> x*y*z,
        (x, y, z) -> x^2*y^2*z^2,
        (x, y, z) -> log(A*x^2 + B*y^2 + C*z^2 + one(x)) 
    ]

    funs_x = [
        (x, y, z) -> one(x),
        (x, y, z) -> y*z,
        (x, y, z) -> 2x*y^2*z^2,
        (x, y, z) -> 2*A*x/(A*x^2 + B*y^2 + C*z^2 + one(x))
    ]

    funs_y = [
        (x, y, z) -> one(x),
        (x, y, z) -> x*z,
        (x, y, z) -> x^2*2y*z^2,
        (x, y, z) -> 2*B*y/(A*x^2 + B*y^2 + C*z^2 + one(x))
    ]

    funs_z = [
        (x, y, z) -> one(x),
        (x, y, z) -> x*y,
        (x, y, z) -> x^2*y^2*2z,
        (x, y, z) -> 2*C*z/(A*x^2 + B*y^2 + C*z^2 + one(x))
    ]

    funs_xy = [
        (x, y, z) -> zero(x),
        (x, y, z) -> z,
        (x, y, z) -> 4*x*y*z^2,
        (x, y, z) -> -4*A*B*x*y/(A*x^2 + B*y^2 + C*z^2 + one(x))^2
    ]

    funs_xz = [
        (x, y, z) -> zero(x),
        (x, y, z) -> y,
        (x, y, z) -> 4*x*y^2*z,
        (x, y, z) -> -4*A*C*x*z/(A*x^2 + B*y^2 + C*z^2 + one(x))^2
    ]

    funs_yz = [
        (x, y, z) -> zero(x),
        (x, y, z) -> x,
        (x, y, z) -> 4*x^2*y*z,
        (x, y, z) -> -4*B*C*y*z/(A*x^2 + B*y^2 + C*z^2 + one(x))^2
    ]

    funs_xyz = [
        (x, y, z) -> zero(x),
        (x, y, z) -> one(x),
        (x, y, z) -> 8*x*y*z,
        (x, y, z) -> 16*A*B*C*x*y*z/(A*x^2 + B*y^2 + C*z^2 + one(x))^3
    ]

    @testset "3D Basic" begin
        @out "\n3D Basic"
        @testset "$T" for T in [Float32 Float64]
        @out "$T" 1
            #Test interval

            xs  = 1:T(0.2):2
            xs0 = xs[1:end-1]
            xst = xs0 + step(xs)/3

            ys  = 5:T(0.4):7
            ys0 = ys[1:end-1]'
            yst = ys0 + step(ys)/3

            zs  = 2:T(1):7
            zs0 = Array{T}(1, 1, length(zs) - 1)
            zs0[1, 1, :] = zs[1:end-1]
            zst = zs0 + step(zs)/3
            
            @testset "$b" for b in instances(Boundary)
            @out "$b" 2

                @testset "$f" for (f, f_x, f_y, f_z, f_xy, f_xz, f_yz, f_xyz) in 
                            zip(funs, funs_x, funs_y, funs_z, 
                                funs_xy, funs_xz, funs_yz, funs_xyz)
                    @out "$f" 3

                    ip          = Interpolation(xs, ys, zs, f, f_x, f_y,
                                  f_z, f_xy,f_xz, f_yz, f_xyz, b)

                    ipf(x,y,z)    = interpolate(ip, x, y, z)

                    ip_x(x,y,z)   = diff_x(ip, x, y, z)
                    ip_y(x,y,z)   = diff_y(ip, x, y, z)
                    ip_z(x,y,z)   = diff_z(ip, x, y, z)

                    ip_xy(x,y,z)  = diff(ip, x, y, z, 1, 1, 0)
                    ip_xz(x,y,z)  = diff(ip, x, y, z, 1, 0, 1)
                    ip_yz(x,y,z)  = diff(ip, x, y, z, 0, 1, 1)
                    ip_xyz(x,y,z) = diff(ip, x, y, z, 1, 1, 1)

                    @test isapprox(ipf.(xst, yst, zst),ip.(xst, yst, zst))

                    @testset "interpolation at xp = 0" begin
                        @test isapprox(ip.(xs0, ys0, zs0), 
                                       f.(xs0, ys0, zs0))
                        @test isapprox(ip_x.(xs0, ys0, zs0), 
                                       f_x.(xs0, ys0, zs0))
                        @test isapprox(ip_y.(xs0, ys0, zs0), 
                                       f_y.(xs0, ys0, zs0))
                        @test isapprox(ip_z.(xs0, ys0, zs0), 
                                       f_z.(xs0, ys0, zs0))
                        @test isapprox(ip_xy.(xs0, ys0, zs0), 
                                       f_xy.(xs0, ys0, zs0), atol=tols[3]) 
                        @test isapprox(ip_xz.(xs0, ys0, zs0),   
                                       f_xz.(xs0, ys0, zs0), atol=tols[3]) 
                        @test isapprox(ip_yz.(xs0, ys0, zs0), 
                                       f_yz.(xs0, ys0, zs0), atol=tols[3]) 
                        @test isapprox(ip_xyz.(xs0, ys0, zs0),
                                       f_xyz.(xs0, ys0, zs0), atol=tols[4]) 
                    end

                    @testset "interpolation in cells" begin
                        @test isapprox(ip.(xst, yst, zst), f.(xst, yst, zst),
                                       atol=tols[2], rtol=tols[1])
                        @test isapprox(ip_x.(xst, yst, zst), 
                                       f_x.(xst, yst, zst),
                                       atol=tols[3], rtol=tols[2])
                        @test isapprox(ip_y.(xst, yst, zst), 
                                       f_y.(xst, yst, zst),
                                       atol=tols[3], rtol=tols[2])
                        @test isapprox(ip_z.(xst, yst, zst), 
                                       f_z.(xst, yst, zst),
                                       atol=tols[3], rtol=tols[2])
                        @test isapprox(ip_xy.(xst, yst, zst), 
                                       f_xy.(xst, yst, zst),
                                       atol=tols[4], rtol=tols[3])
                        @test isapprox(ip_xz.(xst, yst, zst), 
                                       f_xz.(xst, yst, zst),
                                       atol=tols[4], rtol=tols[3])
                        @test isapprox(ip_yz.(xst, yst, zst), 
                                       f_yz.(xst, yst, zst),
                                       atol=tols[4], rtol=tols[3])
                        @test isapprox(ip_xyz.(xst, yst, zst),
                                       f_xyz.(xst, yst, zst),
                                       atol=tols[4], rtol=tols[4])
                    end
                end
            end
        end
    end
end

function exp3d()
    tols= [1e-3, 1e-2, 1e-1, 1]

    xs  = 1:0.2:2
    xs0 = xs[1:end-1]
    xst = xs0 + step(xs)/3

    ys  = 5:0.4:7
    ys0 = ys[1:end-1]'
    yst = ys0 + step(ys)/3

    zs  = 2:1.:7
    zs0 = Array{Float64}(1, 1, length(zs) - 1)
    zs0[1, 1, :] = zs[1:end-1]
    zst = zs0 + step(zs)/3

    exp3d(x, y, z) = exp(x + y + z)

    ip  = Interpolation(xs, ys, zs, ntuple(i->exp3d, 8)...)

    @testset "3D exp" begin
        @out "\n3D exp"
        @test isapprox(ip.(xs0, ys0, zs0), exp3d.(xs0, ys0, zs0))
        @test isapprox(ip.(xst, yst, zst), exp3d.(xst, yst, zst), 
                       rtol=tols[2], atol=tols[2])

        for i = 0:2
            for j = 0:2
                for k = 0:2
                    n   = i+j+k
                    if n > 0
                        m   = maximum([i, j, k]) + 1
                        @out "d^$(n)f/dx^$i/dy^$j/dy^$k" 1
                        dip(x, y, z)  = diff(ip, x, y, z, i, j, k)
                        @test isapprox(dip.(xst, yst, zst), 
                                       exp3d.(xst, yst, zst), 
                                       rtol=tols[m], atol=tols[m])
                    end
                end
            end
        end


    end
end

function boundary3d()
    tols= [1e-3, 1e-2, 1e-1, 1]

    fr  = 89//97
    big = 1e6

    xs = 0:0.1:1
    ys = 5:0.2:7
    zs = 2:0.5:7

    A   = 2
    B   = 1//2
    C   = 3//4
    f(x, y, z)      = log(A*x^2 + B*y^2 + C*z^2 + one(x)) 
    f_x(x, y, z)    = 2*A*x/(A*x^2 + B*y^2 + C*z^2 + one(x))
    f_y(x, y, z)    = 2*B*x/(A*x^2 + B*y^2 + C*z^2 + one(x))
    f_z(x, y, z)    = 2*C*x/(A*x^2 + B*y^2 + C*z^2 + one(x))
    f_xy(x, y, z)   = -4*A*B*x*y/(A*x^2 + B*y^2 + C*z^2 + one(x))^2
    f_xz(x, y, z)   = -4*A*C*x*y/(A*x^2 + B*y^2 + C*z^2 + one(x))^2
    f_yz(x, y, z)   = -4*B*C*x*y/(A*x^2 + B*y^2 + C*z^2 + one(x))^2
    f_xyz(x, y, z)  = 16*A*B*C*x*y*z/(A*x^2 + B*y^2 + C*z^2 + one(x))^3

    dx  = xs[end] - xs[1]
    dy  = ys[end] - ys[1]
    dz  = zs[end] - zs[1]

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

    @testset "3D Boundary" begin
        @out "\n3D Boundary"
        @testset "Mode throw_error" begin
            @out "throw_error" 1
           ip = Interpolation(xs, ys, zs, f, f_x, f_y,
                                  f_z, f_xy,f_xz, f_yz, f_xyz, throw_error)

            ip_x(x, y, z)  = diff_x(ip, x, y, z)
            ip_y(x, y, z)  = diff_y(ip, x, y, z)
            ip_z(x, y, z)  = diff_z(ip, x, y, z)

            ip_xx(x, y, z) = diff(ip, x, y, z, 2, 0, 0)
            ip_xy(x, y, z) = diff(ip, x, y, z, 1, 1, 0)
            ip_xz(x, y, z) = diff(ip, x, y, z, 1, 0, 1)
            ip_yy(x, y, z) = diff(ip, x, y, z, 0, 2, 0)
            ip_yz(x, y, z) = diff(ip, x, y, z, 0, 1, 1)
            ip_zz(x, y, z) = diff(ip, x, y, z, 0, 0, 2)

            for (xt, yt, zt) in collect(Base.product(xts, yts, zts))[2:end]  
                @test_throws BoundsError ip(xt, yt, zt)
                @test_throws BoundsError ip_x(xt, yt, zt)
                @test_throws BoundsError ip_y(xt, yt, zt)
                @test_throws BoundsError ip_z(xt, yt, zt)
                @test_throws BoundsError ip_xx(xt, yt, zt)
                @test_throws BoundsError ip_xy(xt, yt, zt)
                @test_throws BoundsError ip_xz(xt, yt, zt)
                @test_throws BoundsError ip_yy(xt, yt, zt)
                @test_throws BoundsError ip_yz(xt, yt, zt)
                @test_throws BoundsError ip_zz(xt, yt, zt)
            end
        end

        #NO TESTING FOR UNSAFE!

        @testset "Mode undef" begin
            @out "undef" 1

           ip = Interpolation(xs, ys, zs, f, f_x, f_y,
                                  f_z, f_xy,f_xz, f_yz, f_xyz, undef)

            ip_x(x, y, z)  = diff_x(ip, x, y, z)
            ip_y(x, y, z)  = diff_y(ip, x, y, z)
            ip_z(x, y, z)  = diff_z(ip, x, y, z)

            ip_xx(x, y, z) = diff(ip, x, y, z, 2, 0, 0)
            ip_xy(x, y, z) = diff(ip, x, y, z, 1, 1, 0)
            ip_xz(x, y, z) = diff(ip, x, y, z, 1, 0, 1)
            ip_yy(x, y, z) = diff(ip, x, y, z, 0, 2, 0)
            ip_yz(x, y, z) = diff(ip, x, y, z, 0, 1, 1)
            ip_zz(x, y, z) = diff(ip, x, y, z, 0, 0, 2)
            
            for (xt, yt, zt) in collect(Base.product(xts, yts, zts))[2:end]  
                @test isnan(ip(xt, yt, zt))
                @test isnan(ip_x(xt, yt, zt))
                @test isnan(ip_y(xt, yt, zt))
                @test isnan(ip_z(xt, yt, zt))
                @test isnan(ip_xx(xt, yt, zt))
                @test isnan(ip_xy(xt, yt, zt))
                @test isnan(ip_xz(xt, yt, zt))
                @test isnan(ip_yy(xt, yt, zt))
                @test isnan(ip_yz(xt, yt, zt))
                @test isnan(ip_zz(xt, yt, zt))
            end
        end
    end
end
