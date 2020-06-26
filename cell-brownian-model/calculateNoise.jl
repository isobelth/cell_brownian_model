#
#  calculateNoise.jl
#  cell-brownian-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module CalculateNoise

using Random
using Distributions
using LinearAlgebra
using StaticArrays

@inline function calculateNoise!(W::MMatrix,Ncells::Int64,dt::Float64)

    # Loop over all trimers
    for ii=1:Ncells
        #ranTheta = 2.0*pi*rand(Uniform(0.0,1.0))
        ranTheta1 = 2.0 * pi * rand(Uniform(0.0, 1.0))
        ranTheta2 = 2.0 * pi * rand(Uniform(0.0, 1.0))
        ranR = abs(rand(Normal(0.0,sqrt(dt))))

        #W[ii,1] += ranR*cos(ranTheta)
        #W[ii,2] += ranR*sin(ranTheta)
        W[ii,1] += ranR*sin(ranTheta1)*cos(ranTheta2)
        W[ii,2] += ranR*sin(ranTheta1)*sin(ranTheta2)
        W[ii,3] += ranR * cos(ranTheta1)

    end
end

export calculateNoise!

end
