__precompile__()
module InterCellForces

using LinearAlgebra
using StaticArrays
using PyCall
sp = pyimport("scipy.spatial")

@inline function interCellForces!(pos,F,Ncells,ϵ,σ,r,a, age, lifetime, m, vertex_points, CoM, unit_vecs)
#finding the Morse Force between the cell pairs
    for ii=1:Ncells
        for jj=1:Ncells
            if ii==jj
                #skip
            else
                r .= pos[jj,:] .- pos[ii,:]
                r_mag = sqrt(dot(r,r))
                #r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-σ)/σ)-exp(-2.0*a*(r_mag-σ)/σ))/r_mag;
                equilibrium_sep =σ*(1+ (0.13*(age[jj]+age[ii]))/(2*lifetime))
                r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-equilibrium_sep))-exp(-2.0*a*(r_mag-equilibrium_sep)))/r_mag;
                #r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-σ))-exp(-2.0*a*(r_mag-σ)))/r_mag;
                F[ii,:] .+= r
                F[jj,:] .-= r
            end
        end
    end

#linear force for cells F=kx, where x is the distance between the cell and the CoM
    for i in vertex_points
        distance = pos[i,:] .- CoM
        distance_mag = sqrt(dot(distance,distance))
        F[i,:] .+= (unit_vecs[i,:] .* m)
    end

end

export interCellForces!

end
