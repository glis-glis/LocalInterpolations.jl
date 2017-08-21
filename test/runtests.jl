################################################################################
#  Copyright (C) 2017 Andreas Füglistaler <andreas.fueglistaler@gmail.com>
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
################################################################################


using LocalInterpolations
using Base.Test

verbalize = length(ARGS) > 0 && ARGS[1] == "-v"

macro out(s, l=0)
    if verbalize
        return :(println(repeat(" ", 2*$l), $(esc(s))))
    else
        return :()
    end
end

include("1D.jl")
include("2D.jl")
include("3D.jl")
include("4D.jl")

tic()
@out "1D"
basic1d()
exp1d()
boundary1d()

@out "2D"
basic2d()
exp2d()
boundary2d()

@out "3D"
basic3d()
exp3d()
boundary3d()

@out "4D"
basic4d()
exp4d()
boundary4d()

t   = toq()
@out "Elapsed time: $(t)s"
