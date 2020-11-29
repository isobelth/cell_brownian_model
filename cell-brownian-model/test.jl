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
    Ncells=10
    age            = MMatrix{Ncells,1}(zeros(Ncells,1))
    pos            = MMatrix{Ncells,3}(zeros(Ncells,3))
    posBM = MMatrix{Ncells,3}(zeros(Ncells,3))


    for ii=1:Ncells
        pos[ii,:] .= rand(Uniform(-2,2),3)
        posBM[ii,:] .= rand(Uniform(-2,2),3)
    end

    lifetime::Float64 = 5.0


    #hull = sp.ConvexHull(reshape(filter(!iszero, pos), (Ncells,3)))




    hull = sp.Delaunay(reshape(filter(!iszero, pos), (Ncells,3)))
    #pybuiltin(:isinstance)
    if pybuiltin(:isinstance)(hull, sp.Delaunay) == 0
        hull = sp.Delaunay(hull)
    end
    result = hull.find_simplex(reshape(filter(!iszero, posBM), (Ncells,3)))





    #return neighbours, pos, unit_vecs
    return pos, posBM, result, result.>= 0

end


print(main_test())
#print(A)
