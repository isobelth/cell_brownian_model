# Import Julia packages
using Distributions
using StaticArrays
using PyCall
using Revise
sp = pyimport("scipy.spatial")

# Import program modules
push!(LOAD_PATH,"./")
using InterCellForces
using CalculateNoise
using UpdateSystem
using OutputData
using Initialise
using CreateRunDirectory
using ImportantParameters
using BMForces


# Define function for bringing together modules to run simulation
function main()
    # Define run parameters
    Ncells        = 8           # Number of cells to start with
    σ              = 0.5         # Minimum diameter of cell/Equilibrium separation
    boxSize        = 3.0         # Dimensions of cube in which particles are initialised
    Ncells_max     = 20          # max number of cells

    N_BM_cells = 5
    N_BM_cells_max = 100
    σ_BM              = 0.2        # BM diameter

    # Thermodynamic parameters
    μ              = 1.0           # Fluid viscosity
    kT             = 0.5           # Boltzmann constant*Temperature

    # Force parameters for cells
    ϵ              = 10.0*kT       # Potential depth
    k              = 10.0*kT       # Cell stiffness/approximate equilibrium spring constant of morse potential
    m             = 5.0            # proportionality constant for matrix force

    # Derived parameters
    D              = kT/(6.0*π*μ*σ)# Diffusion constant
    a              = sqrt(k/(2*ϵ)) # Property of morse potential derived from approximate equilibrium spring constant

    # Simulation parameters
    tmax           = 10.0          # Total simulation time
    dt             = 0.0001        # Time step between iterations
    outputInterval = tmax/100.0    # Time interval for writing position data to file

    lifetime = 8.0

    # Data arrays
    pos            = MMatrix{Ncells_max,3}(zeros(Ncells_max,3)) # xyz positions of all particles
    F              = MMatrix{Ncells_max,3}(zeros(Ncells_max,3)) # xyz dimensions of all forces applied to particles
    W              = MMatrix{Ncells_max,3}(zeros(Ncells_max,3)) # xyz values of stochastic Wiener process for all particles
    age            = MMatrix{Ncells_max,1}(zeros(Ncells_max,1))
    BM_pos         = MMatrix{N_BM_cells_max,3}(zeros(N_BM_cells_max,3))
    FBM            = MMatrix{N_BM_cells_max,3}(zeros(N_BM_cells_max,3))
    WBM            = MMatrix{N_BM_cells_max,3}(zeros(N_BM_cells_max,3))

    # Allocate variables needed for calculations
    t = 0.0
    AA = zeros(3)
    BB = zeros(3)
    CC  = zeros(4)
    DD = MArray{Tuple{3},Float64,1,3}(zeros(3))

    # Set up folder and output files for data
    foldername = createRunDirectory(Ncells,Ncells_max, lifetime,boxSize, k,μ,kT,ϵ,σ,D,tmax,dt,outputInterval, m)
    outfile = open("output/$(foldername)/celloutput.txt","w")
    outfile2 = open("output/$(foldername)/BM.txt","w")


    # Initialise cell locations within box
    initialise(pos,Ncells,boxSize, age, lifetime, N_BM_cells, BM_pos)

    # Iterate over system time
    while t<tmax

        # Save position data to file
        if (t%outputInterval)<dt
            outputData(pos,outfile,outfile2, t,tmax, age, BM_pos)
        end

        #create memory for the unit vecs, COM and neighbours
        unit_vecs=zeros(Ncells, 3)
        neighbours = zeros(Ncells, Ncells)
        CoM::Float64 = 0.0

        unit_vecs_BM=zeros(N_BM_cells, 3)
        neighbours_BM = zeros(N_BM_cells, N_BM_cells)
        CoM_total::Float64 = 0.0

        #calculate the convex hull and the vertex points
        hull = sp.ConvexHull(reshape(filter(!iszero, pos), (Ncells,3)))
        vertex_points = hull.vertices.+1

        print("BM pos = ", BM_pos, "\n")
        hull_BM = sp.ConvexHull(reshape(filter(!iszero, BM_pos), (N_BM_cells,3)))
        vertex_points_BM = hull_BM.vertices.+1

        #function calculates the neighbours list and CoM
        importantParameters!(pos, BM_pos, Ncells, N_BM_cells, age, t, lifetime, unit_vecs, hull, neighbours, CoM, vertex_points, unit_vecs_BM, hull_BM, neighbours_BM, CoM_total, vertex_points_BM)
        #calculates forces between cells
        interCellForces!(pos,F,Ncells,ϵ,σ,DD,a, age, lifetime, m, vertex_points, CoM, unit_vecs)
        #calculates forces between BM cells
        bMForces!(pos,BM_pos,FBM,Ncells, N_BM_cells,ϵ,σ,σ_BM, DD,a,t, lifetime, vertex_points_BM, CoM_total, unit_vecs_BM)

        # Adapt timestep to maximum force value
        Fmax_sq = max(sum(F.*F,dims=2)...) #Fdims from 2
        dt = min(σ^2/(4.0*32*D),kT*σ/(4.0*D*sqrt(Fmax_sq)))

        # Calculate stochastic component
        calculateNoise!(W,Ncells,WBM, N_BM_cells,dt)

        # Forward euler integration of system
        t, Ncells, N_BM_cells = updateSystem!(DD, pos,F,W,Ncells, BM_pos, FBM, WBM, N_BM_cells, t,dt,D,kT,age,lifetime,σ,σ_BM, hull, vertex_points, hull_BM, vertex_points_BM, neighbours_BM)
        print("t=",t,"Ncells=",Ncells,"NBMcells=",N_BM_cells,"\n")
        #print("NeighboursBM=",neighbours_BM,"\n")

    end
    close(outfile)
end

main()
