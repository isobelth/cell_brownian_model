using StaticArrays
@pyimport matplotlib.pyplot as pyplt
using PyCall
#using VoronoiDelaunay
sp = pyimport("scipy.spatial")
using Distributions
using Plots
using PyPlot
using LinearAlgebra

function main_test()
    Ncells=10
    age            = MMatrix{Ncells,1}(zeros(Ncells,1))
    pos            = MMatrix{Ncells,3}(zeros(Ncells,3))
    for ii=1:Ncells
        pos[ii,:] .= rand(Uniform(-2,2),3)
    end

    lifetime::Float64 = 5.0

    numerator = 0
    denominator::Float64 = 0

    for i in 1:Ncells
        factor::Float64 = (1+ (age[i]/lifetime))^3
        numerator += pos[i,:] .* factor
        denominator += factor
    end

    A = numerator/denominator

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
        print("i=", i, "The neighbours of i are",neighbours[i,:], "\n")
        for j in neighbours[i,:]
            if j != 0.0
                j=convert(Int64, j)
                vector2 = pos[j,:] .- pos[i,:]
                mag_2 = sqrt(dot(vector2, vector2))
                #mag_2 = sqrt(vector2.dot(vector2))
                unit_vector2 = vector2/mag_2
                unit_vector += unit_vector2
            end
        end

        if sqrt(dot(unit_vector, unit_vector)) != 0
            unit_vector = unit_vector / (sqrt(dot(unit_vector,unit_vector)))
            #unit_vector = unit_vector.tolist()
        end

        unit_vecs[i,:] = unit_vector
    end
    print("The unit vectors are", unit_vecs, "\n")
    print("positions are = ", pos, "\n")
    print("neighbours = ", neighbours, "\n")
    proportionality_constant = 1.0
    for i in vertex_points
        distance = pos[i,:].- CoM


    #return neighbours, pos, unit_vecs


end


print(main_test())
#print(A)
