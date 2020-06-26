#
#  interCellForces.jl
#  cell-brownian-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module InterCellForces

using LinearAlgebra
using StaticArrays

@inline function interCellForces!(pos,F,Ncells,ϵ,σ,r,a, age, lifetime,)

    for ii=1:Ncells
        for jj=1:Ncells
            if ii==jj
                #skip
            else
                r .= pos[jj,:] .- pos[ii,:]
                r_mag = sqrt(dot(r,r))
                #r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-σ)/σ)-exp(-2.0*a*(r_mag-σ)/σ))/r_mag;
                equilibrium_sep =0.5*σ*(1+ (age[jj]+age[ii])/(2*lifetime))
                r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-equilibrium_sep))-exp(-2.0*a*(r_mag-equilibrium_sep)))/r_mag;
                #r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-σ))-exp(-2.0*a*(r_mag-σ)))/r_mag;

                F[ii,:] .+= r
                F[jj,:] .-= r
            end
        end
    end

end

export interCellForces!

end
