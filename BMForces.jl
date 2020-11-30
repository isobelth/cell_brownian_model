#

module BMForces

using LinearAlgebra
using StaticArrays
using PyCall
sp = pyimport("scipy.spatial")

@inline function bMForces!(pos,BM_pos, FBM,Ncells, N_BM_cells,ϵ,σ,σ_BM, r,a,t,lifetime, vertex_points_BM, CoM_total, unit_vecs_BM)
#############MORSE FORCE BETWEEN BMCELLS##########
    for ii=1:N_BM_cells
        for jj=1:N_BM_cells
            if ii==jj
                #skip
            else
                r .= BM_pos[jj,:] .- BM_pos[ii,:]
                r_mag = sqrt(dot(r,r))
                equilibrium_sep =(σ_BM)
                r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-equilibrium_sep))-exp(-2.0*a*(r_mag-equilibrium_sep)))/r_mag;
                #r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-σ))-exp(-2.0*a*(r_mag-σ)))/r_mag;
                FBM[ii,:] .+= r
                FBM[jj,:] .-= r
            end
        end
    end
#########IF THE BM CELL IS INSIDE THE CONVEX HULL, WE ADD LOADS OF FORCE IN THE OUTWARD DIRECTION

    # hull = sp.Delaunay(reshape(filter(!iszero, pos), (Ncells,3)))
    # #pybuiltin(:isinstance)
    # if pybuiltin(:isinstance)(hull, sp.Delaunay) == 0
    #     hull = sp.Delaunay(hull)
    # end
    # result = hull.find_simplex(reshape(filter(!iszero, BM_pos), (N_BM_cells,3)))
    #go through each cell and add 100 (outwards) to its force if its within the BM
    # for i in N_BM_cells
    #     if result[i] != -1
    #         FBM[i,:] .+= (unit_vecs_BM[i,:] .* result .* 100)
    #     end
    # end




    #return neighbours, pos, unit_vecs
    #return pos, posBM, result, result.>= 0

end

export bMForces!

end
