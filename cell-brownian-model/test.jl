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


    #hull = sp.ConvexHull(reshape(filter(!iszero, pos), (Ncells,3)))


    def in_hull(p,hull):
        if not isinstance(hull,Delaunay):
            hull = Delaunay(hull)

            return hull.find_simplex(p)>=0

    bm_cells = np.random.rand(10,3)
    hull = sp.Delaunay(reshape(filter(!iszero, pos), (Ncells,3)))
    result = in_hull(bm_cells, hull)


    #return neighbours, pos, unit_vecs


end


print(main_test())
#print(A)
