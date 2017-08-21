################################################################################
#  Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

#Test interval
function basic2d()
    tols= [1e-3, 1e-2, 1e-1]

    #Test functions
    A   = 2
    B   = 1//2

    funs = [ 
        (x, y) -> x + y,
        (x, y) -> x*y,
        (x, y) -> x^2*y^2,
        (x, y) -> log(A*x^2 + B*y^2 + 1) 
    ]

    funs_x = [
        (x, y) -> one(x),
        (x, y) -> y,
        (x, y) -> 2x*y^2,
        (x, y) -> 2*A*x/(A*x^2 + B*y^2 + 1)
    ]

    funs_y = [
        (x, y) -> one(x),
        (x, y) -> x,
        (x, y) -> 2y*x^2,
        (x, y) -> 2*B*y/(A*x^2 + B*y^2 + 1)
    ]

    funs_xy = [
        (x, y) -> zero(x),
        (x, y) -> one(x),
        (x, y) -> 4*x*y,
        (x, y) -> -4*A*B*x*y/(A*x^2 + B*y^2 + 1)^2
    ]

     

    @testset "2D Basic" begin
        @out "\n2D Basic"
        @testset "$T" for T in [Float32 Float64]
        @out "$T" 1
            xs  = 0:T(0.3):3
            xs0 = xs[1:end-1]
            xst = xs0 + step(xs)/3

            ys  = 5:T(0.2):15
            ys0 = ys[1:end-1]
            yst = ys0 + step(ys)/3

            @testset "$b" for b in instances(Boundary)
            @out "$b" 2

                @testset "$f" for (f, f_x, f_y, f_xy) in 
                            zip(funs, funs_x, funs_y, funs_xy)
                    @out "$f" 3

                    ip          = Interpolation(xs, ys, f, f_x,
                                  f_y, f_xy, b)

                    ipf(x,y)    = interpolate(ip, x, y)
                    ip_x(x,y)   = diff_x(ip, x, y)
                    ip_y(x,y)   = diff_y(ip, x, y)
                    ip_xy(x,y)  = diff(ip, x, y, 1, 1)

                    @test       isapprox(ipf.(xst, yst'),ip.(xst, yst'))

                    @testset "interpolation at xp = 0" begin
                        @test isapprox(ip.(xs0, ys0'), f.(xs0, ys0'))
                        @test isapprox(ip_x.(xs0, ys0'),f_x.(xs0, ys0'))
                        @test isapprox(ip_y.(xs0, ys0'),f_y.(xs0, ys0'))
                        @test isapprox(ip_xy.(xs0, ys0'),f_xy.(xs0, ys0'),
                        atol=tols[3]) # if f_xy = 0
                    end


                    @testset "interpolation in cells" begin
                        @test isapprox(ip.(xst, yst'), f.(xst, yst'), 
                                atol=tols[1], rtol=tols[1])
                        @test isapprox(ip_x.(xst, yst'),f_x.(xst, yst'),
                                atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_y.(xst, yst'),f_y.(xst, yst'),
                                atol=tols[2], rtol=tols[2])
                        @test isapprox(ip_xy.(xst, yst'),f_xy.(xst, yst'),
                                atol=tols[3], rtol=tols[3])
                    end
                end
            end
        end
    end
end

function exp2d()
    tols= [1e-3, 1e-2, 1e-1]

    xs  = 0:0.3:3
    xs0 = xs[1:end-1]
    xst = xs0 + step(xs)/3

    ys  = 5:0.2:15
    ys0 = ys[1:end-1]
    yst = ys0 + step(ys)/3

    exp2d(x, y) = exp(x + y)

    ip  = Interpolation(xs, ys, ntuple(i->exp2d, 4)...)

    @testset "2D exp" begin
        @out "\n2D exp"
        @test isapprox(ip.(xs0, ys0'), exp2d.(xs0, ys0'))
        @test isapprox(ip.(xst, yst'), exp2d.(xst, yst'), rtol=tols[1])

        for i = 0:2
            for j = 0:2
                n   = i+j
                if n > 0
                    m   = maximum([i, j]) + 1
                    @out "d^$(n)f/dx^$i/dy^$j" 1
                    dip(x, y)     = diff(ip, x, y, i, j)
                    @test isapprox(dip.(xst, yst'), exp2d.(xst, yst'), 
                            rtol=tols[m])
                end
            end
        end
    end
end

function boundary2d()
    tols= [1e-3, 1e-2, 1e-1]

    fr  = 89//97
    big = 1e6

    xs  = 0:0.3:3
    ys  = 5:0.2:15

    A          = 2
    B          = 1//2
    f(x, y)    = log(A*x^2 + B*y^2 + 1) 
    f_x(x, y)  = 2*A*x/(A*x^2 + B*y^2 + 1)
    f_y(x, y)  = 2*B*x/(A*x^2 + B*y^2 + 1)
    f_xy(x, y) = -4*A*B*x*y/(A*x^2 + B*y^2 + 1)^2

    dx  = xs[end] - xs[1]
    dy  = ys[end] - ys[1]

    xts = [
        xs[1] + dx*fr,      #inbounds
        xs[end]  + 1/big,
        xs[1]    - 1/big,
        xs[end] +  dx*fr,
        xs[1] -  dx*fr,
        xs[end] + big,
        xs[1] - big
    ]

    yts = [
        ys[1] + dy*fr,      #inbounds
        ys[end]  + 1/big,
        ys[1]    - 1/big,
        ys[end] +  dy*fr,
        ys[1] -  dy*fr,
        ys[end] + big,
        ys[1] - big
    ]

    @testset "2D Boundary" begin
        @out "\n2D Boundary"
        @testset "Mode throw_error" begin
        @out "throw_error" 1
            ip          = Interpolation(xs, ys, f, f_x, f_y, f_xy, throw_error)

            for (xt, yt) in collect(Base.product(xts, yts))[2:end]  #skip first
                @test_throws BoundsError ip(xt, yt)
                @test_throws BoundsError diff_x(ip, xt, yt)
                @test_throws BoundsError diff_y(ip, xt, yt)
            end
        end

        #NO TESTING FOR UNSAFE!

        @testset "Mode undef" begin
        @out "undef" 1
            ip          = Interpolation(xs, ys, f, f_x, f_y, f_xy, undef)

            for (xt, yt) in collect(Base.product(xts, yts))[2:end]  #skip first
                @test isnan(ip(xt, yt))
                @test isnan(diff_x(ip, xt, yt))
                @test isnan(diff_y(ip, xt, yt))
            end
        end
    end
end
