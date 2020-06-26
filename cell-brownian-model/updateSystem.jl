#
#  updateSystem.jl
#  cell-brownian-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module UpdateSystem
using Random
using Distributions
using LinearAlgebra
using StaticArrays
include("./outerPoints.jl")
using .OuterPoints

@inline function updateSystem!(pos::MMatrix,F::MMatrix,W::MMatrix,Ncells::Int64,t::Float64,dt::Float64,D::Float64,kT::Float64, age::MMatrix, lifetime::Float64, σ::Float64 )

    pos .+= F.*(dt*D/kT) .+ W.*sqrt(2.0*D)
    F .= 0
    W .= 0

    vertex_indices = outerPoints!(pos)
    print("Vertex indices =", vertex_indices, "\n")
    for ii in 1:Ncells
        if ii in vertex_indices
        else
            age[ii] += dt
        end
    end

    for ii=1: Ncells
        if age[ii] >= lifetime
            Theta = pi * rand(Uniform(0.0, 1.0))
            Phi = 2.0 * pi * rand(Uniform(0.0, 1.0))

            pos[ii, 1] += σ * 0.25 * sin(Theta)*cos(Phi) *0.5
            pos[ii,2] += σ * 0.25 * sin(Theta)*sin(Phi) * 0.5
            pos[ii,3] += σ * 0.25 * cos(Theta) * 0.5
            age[ii] = 0
            Ncells += 1

            pos[Ncells, 1] =  pos[ii,1]-(0.5 *σ * sin(Theta)*cos(Phi) *0.5)
            pos[Ncells,2] = pos[ii,2]-(0.5* σ * sin(Theta)*sin(Phi) * 0.5)
            pos[Ncells,3] = pos[ii,3]-(0.5 * σ *  cos(Theta) * 0.5)

        end
    end



    t += dt
    return t, Ncells

end

export updateSystem!

end
