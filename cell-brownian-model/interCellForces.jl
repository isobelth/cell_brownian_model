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
using PyCall
sp = pyimport("scipy.spatial")

@inline function interCellForces!(pos,F,Ncells,ϵ,σ,r,a, age, lifetime, m)

    for ii=1:Ncells
        for jj=1:Ncells
            if ii==jj
                #skip
            else
                r .= pos[jj,:] .- pos[ii,:]
                r_mag = sqrt(dot(r,r))
                #r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-σ)/σ)-exp(-2.0*a*(r_mag-σ)/σ))/r_mag;
                equilibrium_sep =σ*(1+ (age[jj]+age[ii])/(2*lifetime))
                r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-equilibrium_sep))-exp(-2.0*a*(r_mag-equilibrium_sep)))/r_mag;
                #r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-σ))-exp(-2.0*a*(r_mag-σ)))/r_mag;
                F[ii,:] .+= r
                F[jj,:] .-= r
            end
        end
    end
    numerator = 0
    denominator::Float64 = 0

    for i in 1:Ncells
        factor::Float64 = (1+ (age[i]/lifetime))^3
        numerator += pos[i,:] .* factor
        denominator += factor
    end

    CoM= numerator/denominator

    hull = sp.ConvexHull(reshape(filter(!iszero, pos), (Ncells,3)))
    vertex_points = hull.vertices.+1

    neighbours = zeros(Ncells, Ncells)
    for i in  1:(size(hull.simplices, 1))

        first_number = hull.simplices[i,1]+1
        second_number = hull.simplices[i,2] + 1
        third_number = hull.simplices[i,3] + 1

        if first_number ∉ neighbours[second_number, :]
            for j in 1:Ncells
                if neighbours[second_number,j] == 0
                    neighbours[second_number, j] = first_number
                    break
                end
            end
        end
        if first_number ∉ neighbours[third_number, :]
            for j in 1:Ncells
                if neighbours[third_number,j] == 0
                    neighbours[third_number, j] = first_number
                    break
                end
            end
        end

        if second_number ∉ neighbours[first_number, :]
            for j in 1:Ncells
                if neighbours[first_number,j] == 0
                    neighbours[first_number, j] = second_number
                    break
                end
            end
        end
        if second_number ∉ neighbours[third_number, :]
            for j in 1:Ncells
                if neighbours[third_number,j] == 0
                    neighbours[third_number, j] = second_number
                    break
                end
            end
        end

        if third_number ∉ neighbours[first_number, :]
            for j in 1:Ncells
                if neighbours[first_number,j] == 0
                    neighbours[first_number, j] = third_number
                    break
                end
            end
        end
        if third_number ∉ neighbours[second_number, :]
            for j in 1:Ncells
                if neighbours[second_number,j] == 0
                    neighbours[second_number, j] = third_number
                    break
                end
            end
        end
    end


    unit_vecs=zeros(Ncells, 3)
    for i in vertex_points
        unit_vector = 0
        for j in neighbours[i,:]
            if j != 0.0
                j=convert(Int64, j)
                vector2 = pos[j,:] .- pos[i,:]
                mag_2 = sqrt(dot(vector2, vector2))
                unit_vector2 = vector2/mag_2
                unit_vector += unit_vector2
            end
        end

        if sqrt(dot(unit_vector, unit_vector)) != 0
            unit_vector = unit_vector / (sqrt(dot(unit_vector,unit_vector)))
        end

        unit_vecs[i,:] = unit_vector

        for i in vertex_points
            distance = pos[i,:] .- CoM
            distance_mag = sqrt(dot(distance,distance))
            F[i,:] .+= (unit_vecs[i,:] .* m)
        end


    end

end

export interCellForces!

end
