#
#  interCellForces.jl
#  cell-brownian-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

__precompile__()
module ImportantParameters

using LinearAlgebra
using StaticArrays
using PyCall
sp = pyimport("scipy.spatial")

@inline function importantParameters!(pos, Ncells, σ, age, lifetime, unit_vecs, hull, neighbours, CoM, vertex_points)
###################FIND THE CENTRE OF MASS################
    numerator = 0
    denominator::Float64 = 0

    for i in 1:Ncells
        factor::Float64 = (1+ 0.26*(age[i]/lifetime))^3
        numerator += pos[i,:] .* factor
        denominator += factor
    end

    CoM= numerator/denominator
###################FIND THE CONVEX HULL AND LIST OF NEIGHBOURS################
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

############FIND AN APPROXIMATE INWARD DIRECTION##########
    #unit_vecs=zeros(Ncells, 3)
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
    end
end

export interCellForces!

end
