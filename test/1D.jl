################################################################################
#  Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################

function basic1d()
    #Test functions
    funs    = [ 
        one,
        x -> x,
        x -> 2x + one(x),
        x -> log(x^2 + 1) 
    ]

    funs_x  =[
        zero,
        one,
        x -> 2*one(x),
        x -> 2x/(x^2 + 1)
    ]

    #Tolerances for interpolation and derivatives
    tols= [1e-4, 1e-3, 1e-2]

    @testset "1D Basic" begin
        @out "\n1D Basic"
        @testset "$T" for T in [Float32 Float64]
        @out "$T" 1
            xs  =  1:T(0.3):7
            xs0 = xs[1:end-1]
            xst = xs0 + step(xs)/3

            @testset "$b" for b in instances(Boundary)
            @out "$b" 2
                @testset "$f" for (f, f_x) in zip(funs, funs_x)
                    @out "$f" 3
                    ip      = Interpolation(xs, f, f_x, b)
                    #ipf(x)  = interpolate(ip, x)    #why does that not work?
                    #@test isapprox(ipf.(xst),ip.(xst))

                    ip_x(x) = diff_x(ip, x)
        
                    @testset "interpolation at xp = 0" begin
                        @test isapprox(ip.(xs0), f.(xs0))
                        @test isapprox(ip_x.(xs0),f_x.(xs0))
                    end
                    @out "xp=0" 3
                    
                    @testset "interpolation in cells" begin
                        @test isapprox(ip.(xst), f.(xst), rtol=tols[1])
                        @test isapprox(ip_x.(xst),f_x.(xst),rtol=tols[2])
                    end
                end
            end
        end
    end
end


function exp1d()
    tols= [1e-4, 1e-3, 1e-2]

    xs  = 1:0.3:7
    xs0 = xs[1:end-1]
    xst = xs0 + step(xs)/3

    ip          = Interpolation(xs, exp, exp)

    @testset "1D exp" begin
        @out "\n1D exp"
        @test isapprox(ip.(xs0), exp.(xs0))
        @test isapprox(ip.(xst), exp.(xst), rtol=tols[1])

        for i = 1:2
            @out "d^$(i)f/dx^$i" 1
            dip(x)     = diff(ip, x, i)
            @test isapprox(dip.(xst), exp.(xst), rtol=tols[i+1])
        end
    end
end
        
function boundary1d()
    tols= [1e-4, 1e-3, 1e-2]

    big = 1e6
    fr  = 89//97

    xs      = 1:0.3:7
    f(x)    = log(x^2 + 1) 
    f_x(x)  = 2x/(x^2 + 1)

    dx  = xs[end] - xs[1]

    xts = [
        xs[end]  + 1/big,
        xs[1]    - 1/big,
        xs[end] +  dx*fr,
        xs[1] -  dx*fr,
        xs[end] + big,
        xs[1] - big
    ]
    @testset "1D Boundary" begin
        @out "\n1D Boundary"
        @testset "Mode throw_error" begin
        @out "throw_error" 1
            ip  = Interpolation(xs, f, f_x, throw_error)

            for xt in xts
                @test_throws BoundsError ip(xt)
                @test_throws BoundsError diff_x(ip, xt)
            end
        end

        #NO TESTING FOR UNSAFE!

        @testset "Mode undef" begin
        @out "undef" 1
            ip  = Interpolation(xs, f, f_x, undef)

            for xt in xts
                @test isnan(ip(xt))
                @test isnan(diff_x(ip, xt))
            end
        end

    end
end
