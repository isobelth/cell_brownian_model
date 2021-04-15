#
#  updateSystem.jl
#  cell-brownian-model
#
#

module UpdateSystem

using Random
using Distributions
using LinearAlgebra
using StaticArrays
#include("./outerPoints.jl")
#using .OuterPoints
using PyCall
sp = pyimport("scipy.spatial")


@inline function updateSystem!(r, pos::MMatrix,F::MMatrix,W::MMatrix,Ncells::Int64, BM_pos::MMatrix, FBM::MMatrix, WBM::MMatrix, N_BM_cells::Int64, t::Float64,dt::Float64,D::Float64,kT::Float64, age::MMatrix, lifetime::Float64, σ::Float64,σ_BM::Float64, hull, vertex_points, hull_BM, vertex_points_BM, neighbours_BM )

    pos .+= F.*(dt*D/kT) .+ W.*sqrt(2.0*D)
    F .= 0
    W .= 0

    display(FBM)
    display(WBM)
    BM_pos .+= FBM.*(dt*D/kT) .+ WBM.*sqrt(2.0*D)
    FBM .= 0
    WBM .= 0

    for ii in 1:Ncells
        if ii in vertex_points
            age[ii] += dt
        end
    end

    for ii=1: Ncells
        if age[ii] >= lifetime
            Theta = pi * rand(Uniform(0.0, 1.0))
            Phi = 2.0 * pi * rand(Uniform(0.0, 1.0))

            pos[ii, 1] += σ * 0.5 * sin(Theta)*cos(Phi)
            pos[ii,2] += σ * 0.5 * sin(Theta)*sin(Phi)
            pos[ii,3] += σ * 0.5 * cos(Theta)
            age[ii] = 0
            Ncells += 1

            pos[Ncells, 1] =  pos[ii,1]-(σ * sin(Theta)*cos(Phi))
            pos[Ncells,2] = pos[ii,2]-(σ * sin(Theta)*sin(Phi))
            pos[Ncells,3] = pos[ii,3]- (σ *  cos(Theta) )
        end
    end

    tested = Tuple{Int64,Int64}[] #these are all the ones we've checked already
    for i in vertex_points_BM
        for j in neighbours_BM[i,:]
            if j != 0.0
                j=convert(Int64, j)
                if (i,j) in tested || (j,i) in tested #|| means or
                    #skip
                else
                    push!(tested,(i,j))
                    r .= BM_pos[j,:] .- BM_pos[i,:]
                    r_mag = sqrt(dot(r,r))
                    if r_mag > 15*σ_BM
                        N_BM_cells += 1
                        #in each timestep, the positions are given and look to be correct
                        #however, in the next timestep, the positions calculated in the previous timestep are forgotten and replaced with NaN
                        BM_pos[N_BM_cells, :] .= 0.5*(BM_pos[j,:] + BM_pos[i,:])
                    else
                    # skip
                    end
                end
            end
        end
    end

    t += dt
    return t, Ncells, N_BM_cells

end

export updateSystem!

end
