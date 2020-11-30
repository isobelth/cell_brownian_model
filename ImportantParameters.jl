#  cell-brownian-model
#


module ImportantParameters

using LinearAlgebra
using StaticArrays
using PyCall
sp = pyimport("scipy.spatial")

@inline function importantParameters!(pos,BM_pos, Ncells,N_BM_cells, age,t, lifetime, unit_vecs, hull, neighbours, CoM, vertex_points, unit_vecs_BM, hull_BM, neighbours_BM, CoM_total, vertex_points_BM)
###################FIND THE CENTRE OF MASS################
    numerator = 0
    denominator::Float64 = 0
    factor::Float64 = 0
    factor_BM::Float64 = 0
    numerator_total = 0
    denominator_total = 0

#TO FIND COM JUST OF THE CELLS
    for i in 1:Ncells
        factor = (1+ 0.26*(age[i]/lifetime))^3
        numerator += pos[i,:] .* factor
        denominator += factor
        numerator_total += pos[i,:] .*factor
        denominator_total += factor
    end
#TO FIND COM OF THE ENTIRE CLUSTER
    for i in 1:N_BM_cells
        factor_BM = (1 + (t/1000))^3############UPDATE THIS!
        numerator_total += BM_pos[i,:] .*factor_BM
        denominator_total += factor_BM
    end


#COM of just the cells
    CoM= numerator/denominator
#COM of the cells and the BM
    CoM_total = numerator_total / denominator_total
###################FIND THE CONVEX HULL AND LIST OF NEIGHBOURS################
###############THIS PART IS FOR THE ACTUAL CELLS#########
    for i in  1:(size(hull.simplices, 1))
        first_number = hull.simplices[i,1] .+1 #######
        second_number = hull.simplices[i,2] .+ 1#######
        third_number = hull.simplices[i,3] .+ 1########

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
###########FIND CONVEX HULL AND NEIGHBOURS FOR BM CELLS######################
############ISSUE - this is now just for the BM cells - what if a normal cell is on the outside?
    for i in  1:(size(hull_BM.simplices, 1))
        first_number = hull_BM.simplices[i,1].+1#########
        second_number = hull_BM.simplices[i,2] .+ 1#####
        third_number = hull_BM.simplices[i,3] .+ 1######

        if first_number ∉ neighbours_BM[second_number, :]
            for j in 1:N_BM_cells
                if neighbours_BM[second_number,j] == 0
                    neighbours_BM[second_number, j] = first_number
                    break
                end
            end
        end
        if first_number ∉ neighbours_BM[third_number, :]
            for j in 1:N_BM_cells
                if neighbours_BM[third_number,j] == 0
                    neighbours_BM[third_number, j] = first_number
                    break
                end
            end
        end

        if second_number ∉ neighbours_BM[first_number, :]
            for j in 1:N_BM_cells
                if neighbours_BM[first_number,j] == 0
                    neighbours_BM[first_number, j] = second_number
                    break
                end
            end
        end
        if second_number ∉ neighbours_BM[third_number, :]
            for j in 1:N_BM_cells
                if neighbours_BM[third_number,j] == 0
                    neighbours_BM[third_number, j] = second_number
                    break
                end
            end
        end

        if third_number ∉ neighbours_BM[first_number, :]
            for j in 1:N_BM_cells
                if neighbours_BM[first_number,j] == 0
                    neighbours_BM[first_number, j] = third_number
                    break
                end
            end
        end
        if third_number ∉ neighbours_BM[second_number, :]
            for j in 1:N_BM_cells
                if neighbours_BM[second_number,j] == 0
                    neighbours_BM[second_number, j] = third_number
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
                unit_vector += unit_vector2 ###########
            end
        end
        if sqrt(dot(unit_vector, unit_vector)) != 0
            unit_vector = unit_vector / (sqrt(dot(unit_vector,unit_vector)))
        end
        unit_vecs[i,:] = unit_vector
    end


#unit_vecs=zeros(Ncells, 3)
# for i in vertex_points_BM
#     unit_vector = 0
#     for j in neighbours_BM[i,:]
#         if j != 0.0
#             j=convert(Int64, j)
#             vector2 = BM_pos[j,:] .- BM_pos[i,:]
#             mag_2 = sqrt(dot(vector2, vector2))
#             unit_vector2 = vector2/mag_2
#             unit_vector += unit_vector2
#         end
#     end
#     if sqrt(dot(unit_vector, unit_vector)) != 0
#         unit_vector = unit_vector / (sqrt(dot(unit_vector,unit_vector)))
#     end
#     unit_vecs_BM[i,:] = unit_vector
# end
end

export importantParameters!

end
