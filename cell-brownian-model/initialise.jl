#
#  initialise.jl
#  cell-brownian-model
#
#  Created by Christopher Revell on 01/05/2020.
#
#
__precompile__()
module Initialise

using StaticArrays
using Random
using Distributions

@inline function initialise(pos::MMatrix,Ncells::Int64,boxsize::Float64,age::MMatrix,lifetime::Float64, N_BM_cells::Int64, BM_pos::MMatrix)

    # for ii=1:Ncells
    #     pos[ii,:] .= rand(Uniform(-boxsize,boxsize),3) #from 2
    #     age[ii] = rand(Uniform(0, lifetime))
    # end
    for ii = 1:Ncells
        pos[ii,:] .= rand(Uniform(-boxsize,boxsize),3) #from 2
        age[ii] = rand(Uniform(0, lifetime))
        BM_pos[ii,:] .= rand(Uniform(-boxsize, boxsize),3) ###might need to make this larger
    end

    for ii = (Ncells + 1):N_BM_cells
        BM_pos[ii,:] .= rand(Uniform(-boxsize, boxsize),3)
    end

end

export initialise

end
