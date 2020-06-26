#
#  createRunDirectory.jl
#  cell-brownian-model
#
#  Created by Christopher Revell on 02/05/2020.
#
#
__precompile__()
module CreateRunDirectory

using Dates

function createRunDirectory(Ncells,Ncells_max, lifetime,boxSize,k,μ,kT,ϵ,σ,D,tmax,dt,outputInterval)
    date = Dates.format(Dates.now(), "dd-mm-yyyy")
    foldername = "$date,N-$Ncells,Nmax-$Ncells_max,tmax-$tmax,kT-$kT"
    #foldername = (Dates.format(Dates.now(),"dd-mm-yyyy"+"N-$Ncells,tmax-$tmax,a-$a,kT-$kT")
    mkdir("output/$(foldername)")
    open("output/$(foldername)/conditions.txt","w") do conditionsfile
        println(conditionsfile,"Ncells         $Ncells      ")
        println(conditionsfile,"Nmax              $Ncells_max    ")
        println(conditionsfile,"Box size       $boxSize      ")
        println(conditionsfile,"k               $k      ")
        println(conditionsfile,"μ              $μ              ")
        println(conditionsfile,"kT             $kT             ")
        println(conditionsfile,"ϵ              $ϵ            ")
        println(conditionsfile,"σ              $σ              ")
        println(conditionsfile,"D              $D              ")
        println(conditionsfile,"tmax           $tmax           ")
        println(conditionsfile,"dt             $dt             ")
        println(conditionsfile,"outputInterval $outputInterval ")
        println(conditionsfile,"lifetime       $lifetime ")
    end

    return foldername
end

export createRunDirectory

end
