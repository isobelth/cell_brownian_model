#
#  interCellForces.jl
#  cell-brownian-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module BMForces

using LinearAlgebra
using StaticArrays
using PyCall
sp = pyimport("scipy.spatial")

@inline function interCellForces!(pos,FBM,N_BM_cells,ϵ,σ,r,a,t,lifetime)
#############MORSE FORCE BETWEEN BMCELLS##########
    for ii=1:N_BM_cells
        for jj=1:N_BM_cells
            if ii==jj
                #skip
            else
                r .= pos[jj,:] .- pos[ii,:]
                r_mag = sqrt(dot(r,r))
                equilibrium_sep =(σ/5)*(1+ (t/1000))
                r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-equilibrium_sep))-exp(-2.0*a*(r_mag-equilibrium_sep)))/r_mag;
                #r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-σ))-exp(-2.0*a*(r_mag-σ)))/r_mag;
                FBM[ii,:] .+= r
                FBM[jj,:] .-= r
            end
        end
    end





    end

end

export bmForces!

end
