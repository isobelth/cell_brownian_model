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

@inline function interCellForces!(pos,F,Ncells,ϵ,σ,r,a, age, lifetime,)

    for ii=1:Ncells
        for jj=1:Ncells
            if ii==jj
                #skip
            else
                r .= pos[jj,:] .- pos[ii,:]
                r_mag = sqrt(dot(r,r))
                #r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-σ)/σ)-exp(-2.0*a*(r_mag-σ)/σ))/r_mag;
                equilibrium_sep =0.5*σ*(1+ (age[jj]+age[ii])/(2*lifetime))
                r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-equilibrium_sep))-exp(-2.0*a*(r_mag-equilibrium_sep)))/r_mag;
                #r .= r.*2.0*ϵ*a*(exp(-a*(r_mag-σ))-exp(-2.0*a*(r_mag-σ)))/r_mag;

                F[ii,:] .+= r
                F[jj,:] .-= r

                #filled_positions = reshape(filter(!iszero, pos), (Ncells,3))
                numerator = 0
                denominator::Float64 = 0

                for i in 1:Ncells
                    factor::Float64 = (1+ (age[i]/lifetime))^3
                    numerator += pos[i,:] .* factor
                    denominator += factor
                end

                print("Numerator =", numerator, "denominator = ", denominator , "\n")
                CoM = numerator/denominator

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




            end
        end
    end

end

export interCellForces!

end
