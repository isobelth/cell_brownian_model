using StaticArrays
@pyimport matplotlib.pyplot as pyplt
#using VoronoiDelaunay
sp = pyimport("scipy.spatial")
using Distributions
using Plots
using PyPlot

function main()
    Ncells=10
    age            = MMatrix{Ncells,1}(zeros(Ncells,1))
    pos            = MMatrix{Ncells,2}(zeros(Ncells,2))
    for ii=1:(Ncells-2)
        pos[ii,:] .= rand(Uniform(-2,2),2)
    end
    my_list = [1,2,5,6]
    if 6 in my_list
        print("yes")
    end
    # A = filter(!iszero, pos)
    # reshape(A, (Ncells-2, 2))
    # reshape(filter(!iszero, pos), (Ncells-2,2))

    #hull = sp.ConvexHull(pos)
    # print("position =", pos, "\n")
    # print("vertices =", (hull.vertices).+1, "\n")
    # print("Vertex points =", pos[hull.vertices.+1,:], "\n")
    #
    # plt.clf()
    # ax=pyplt.gca()
    # ax.set_xlim((-2, 2))
    # ax.set_ylim((-2, 2))
    # ax.set_aspect("equal")
    # plt.plot(pos[:,1], pos[:,2], lw=0, color = :blue, marker="o")
    #
    # plt.plot(pos[hull.vertices.+1,1], pos[hull.vertices.+1,2], lw=0.25, color = :red, marker="o")
    # fig = gcf()
    # display(fig)
    # vertex_points = MArray{Tuple{3},Float64,1,3}(zeros(3))
    # vertex_points = pos[hull.vertices.+1,:]
    # return vertex_points
    return     reshape(filter(!iszero, pos), (Ncells-2,2))
end

print(main())
#print(A)
