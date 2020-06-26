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
include("./outerPoints.jl")
using .InterCellForces
using .CalculateNoise
using .UpdateSystem
using .OutputData
using .Initialise
using .CreateRunDirectory
using .OuterPoints


# Define function for bringing together modules to run simulation
@inline function main()

    # Define run parameters
    Ncells        = 7            # Number of cells to start with
    σ              = 0.5           # Diameter of cell/Equilibrium separation
    boxSize        = 3.0           # Dimensions of cube in which particles are initialised
    Ncells_max     = 10        #max number of cells

    # Thermodynamic parameters
    μ              = 1.0           # Fluid viscosity
    kT             = 1.0           # Boltzmann constant*Temperature

    # Force parameters
    ϵ              = 10.0*kT       # Potential depth
    k              = 10.0*kT       # Cell stiffness/approximate equilibrium spring constant of morse potential

    # Derived parameters
    D              = kT/(6.0*π*μ*σ)# Diffusion constant
    a              = sqrt(k/(2*ϵ)) # Property of morse potential derived from approximate equilibrium spring constant

    # Simulation parameters
    tmax           = 10.0           # Total simulation time
    dt             = 0.0001        # Time step between iterations
    outputInterval = tmax/100.0    # Time interval for writing position data to file

    lifetime = 10.0

    # Data arrays
    pos            = MMatrix{Ncells_max,3}(zeros(Ncells_max,3)) # xyz positions of all particles
    F              = MMatrix{Ncells_max,3}(zeros(Ncells_max,3)) # xyz dimensions of all forces applied to particles
    W              = MMatrix{Ncells_max,3}(zeros(Ncells_max,3)) # xyz values of stochastic Wiener process for all particles
    age            = MMatrix{Ncells_max,1}(zeros(Ncells_max,1))

    # Allocate variables needed for calculations
    t = 0.0
    AA = zeros(3)
    BB = zeros(3)
    CC  = zeros(4)
    DD = MArray{Tuple{3},Float64,1,3}(zeros(3))

    # Set up folder and output files for data
    foldername = createRunDirectory(Ncells,Ncells_max, lifetime,boxSize, k,μ,kT,ϵ,σ,D,tmax,dt,outputInterval)
    outfile = open("output/$(foldername)/output.txt","w")

    # Initialise cell locations within box
    initialise(pos,Ncells,boxSize, age, lifetime)

    # Iterate over system time
    while t<tmax

        # Save position data to file
        if (t%outputInterval)<dt
            outputData(pos,outfile,t,tmax, age)
        end

        # Calculate Morse interactions between cells
        interCellForces!(pos,F,Ncells,ϵ,σ,DD,a, age, lifetime)

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
