using StaticArrays

using PyCall
@pyimport matplotlib.pyplot as pyplt
#using VoronoiDelaunay
sp = pyimport("scipy.spatial")

using Distributions
using Plots
using PyPlot
using LinearAlgebra

function main_test()
    Ncells=5
    Ncells_max = 20
    age            = MMatrix{Ncells,1}(zeros(Ncells,1))
    pos            = MMatrix{Ncells_max,3}(zeros(Ncells_max,3))



    for ii=1:Ncells
        pos[ii,:] .= rand(Uniform(-2,2),3)
    end

    lifetime::Float64 = 5.0
    neighbours = zeros(Ncells, Ncells)
    hull = sp.ConvexHull(reshape(filter(!iszero, pos), (Ncells,3)))
    vertex_points = hull.vertices.+1

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
    vector = 0
    for i in vertex_points
        for j in neighbours[i,:]
            if j != 0.0
                j=convert(Int64, j)
                vector = pos[j,:] .- pos[i,:]
                r_mag = sqrt(vector[1]*vector[1] + vector[2]*vector[2] + vector[3]*vector[3])
                if r_mag > 3
                    Ncells = Ncells + 1
                    pos[N_BM_cells, :] = 0.5*(pos[j,:] + pos[i,:])

                    print("NEW CELL, i=", i, "j=",j, "mag=", r_mag, "POS OF NEW BM CELL=", pos[Ncells,:], "\n")
                else
                    print("NO i=", i, "j=",j,"mag=", r_mag, "\n")
                end
            end
        end
    end


    #return neighbours, pos, unit_vecs
    return pos, Ncells

end


print(main_test())
#print(A)
