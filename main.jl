#
#  main.jl
#  cell-brownian-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

# Import Julia packages
using Distributions
using StaticArrays
using PyCall
sp = pyimport("scipy.spatial")

# Import program modules
include("./outputData.jl")
include("./calculateNoise.jl")
include("./updateSystem.jl")
include("./interCellForces.jl")
include("./initialise.jl")
include("./createRunDirectory.jl")
include("./importantParameters.jl")
#include("./outerPoints.jl")
using .InterCellForces
using .CalculateNoise
using .UpdateSystem
using .OutputData
using .Initialise
using .CreateRunDirectory
using .ImportantParameters
#using .OuterPoints


# Define function for bringing together modules to run simulation
@inline function main()

    # Define run parameters
    Ncells        = 7            # Number of cells to start with
    σ              = 0.5           # Diameter of cell/Equilibrium separation
    boxSize        = 3.0           # Dimensions of cube in which particles are initialised
    Ncells_max     = 70        #max number of cells

    N_BM_cells = 200

    # Thermodynamic parameters
    μ              = 1.0           # Fluid viscosity
    kT             = 0.5           # Boltzmann constant*Temperature

    # Force parameters for cells
    ϵ              = 10.0*kT       # Potential depth
    k              = 10.0*kT       # Cell stiffness/approximate equilibrium spring constant of morse potential
    m             = 5.0           #proportionality constant for matrix force


    # Derived parameters
    D              = kT/(6.0*π*μ*σ)# Diffusion constant
    a              = sqrt(k/(2*ϵ)) # Property of morse potential derived from approximate equilibrium spring constant

    # Simulation parameters
    tmax           = 50.0           # Total simulation time
    dt             = 0.0001        # Time step between iterations
    outputInterval = tmax/100.0    # Time interval for writing position data to file

    lifetime = 8.0

    # Data arrays
    pos            = MMatrix{Ncells_max,3}(zeros(Ncells_max,3)) # xyz positions of all particles
    F              = MMatrix{Ncells_max,3}(zeros(Ncells_max,3)) # xyz dimensions of all forces applied to particles
    W              = MMatrix{Ncells_max,3}(zeros(Ncells_max,3)) # xyz values of stochastic Wiener process for all particles
    age            = MMatrix{Ncells_max,1}(zeros(Ncells_max,1))
    BM_pos         = MMatrix{N_BM_cells,3}(zeros(N_BM_cells,3))
    FBM            = MMatrix{N_BM_cells,3}(zeros(N_BM_cells,3))
    # Allocate variables needed for calculations
    t = 0.0
    AA = zeros(3)
    BB = zeros(3)
    CC  = zeros(4)
    DD = MArray{Tuple{3},Float64,1,3}(zeros(3))

    # Set up folder and output files for data
    foldername = createRunDirectory(Ncells,Ncells_max, lifetime,boxSize, k,μ,kT,ϵ,σ,D,tmax,dt,outputInterval, m)
    outfile = open("output/$(foldername)/output.txt","w")

    # Initialise cell locations within box
    initialise(pos,Ncells,boxSize, age, lifetime, N_BM_cells, BM_pos)

    # Iterate over system time
    while t<tmax

        # Save position data to file
        if (t%outputInterval)<dt
            outputData(pos,outfile,t,tmax, age)
        end

        #calculate CoM, vertices and neighbours list

        # Calculate Morse interactions between cells
        unit_vecs=zeros(Ncells, 3)
        hull = sp.ConvexHull(reshape(filter(!iszero, pos), (Ncells,3)))
        vertex_points = hull.vertices.+1
        neighbours = zeros(Ncells, Ncells)
        CoM::Float64 = 0.0

        #importantParameters!(pos, Ncells, σ, age, lifetime, unit_vecs, hull, neighbours, CoM, vertex_points)

        interCellForces!(pos,F,Ncells,ϵ,σ,DD,a, age, lifetime, m, vertex_points, CoM, unit_vecs)

        bMForces!(pos,FBM,N_BM_cells,ϵ,σ,r,a,t, lifetime)

        # Adapt timestep to maximum force value
        Fmax_sq = max(sum(F.*F,dims=2)...) #Fdims from 2
        dt = min(σ^2/(4.0*32*D),kT*σ/(4.0*D*sqrt(Fmax_sq)))

        # Calculate stochastic component
        calculateNoise!(W,Ncells,dt)

        # Forward euler integration of system
        t, Ncells = updateSystem!(pos,F,W,Ncells,t,dt,D,kT,age,lifetime,σ)
        print("t=",t,"Ncells=",Ncells,"\n")

    end
    close(outfile)
end

main()
