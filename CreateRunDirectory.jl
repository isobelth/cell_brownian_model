#
#  createRunDirectory.jl
#  cell-brownian-model
#
#

module CreateRunDirectory

using Dates

function createRunDirectory(Ncells,Ncells_max, lifetime,boxSize,k,μ,kT,ϵ,σ,D,tmax,dt,outputInterval, m)
    date = Dates.format(Dates.now(), "dd-mm-yyyy")
    foldername = "$date,N-$Ncells,Nmax-$Ncells_max,tmax-$tmax,kT-$kT, m-$m"
    #foldername = (Dates.format(Dates.now(),"dd-mm-yyyy"+"N-$Ncells,tmax-$tmax,a-$a,kT-$kT")
    mkpath("output/$(foldername)")
    open("output/$(foldername)/conditions.txt","w") do conditionsfile
        println(conditionsfile,"Ncells,$Ncells")
        println(conditionsfile,"Nmax,$Ncells_max")
        println(conditionsfile,"Box size,$boxSize")
        println(conditionsfile,"k,$k")
        println(conditionsfile,"μ,$μ")
        println(conditionsfile,"kT,$kT")
        println(conditionsfile,"ϵ,$ϵ")
        println(conditionsfile,"σ,$σ")
        println(conditionsfile,"D,$D")
        println(conditionsfile,"tmax,$tmax")
        println(conditionsfile,"dt,$dt")
        println(conditionsfile,"outputInterval,$outputInterval")
        println(conditionsfile,"lifetime,$lifetime")
        println(conditionsfile,"m,$m")
    end

    return foldername
end

export createRunDirectory

end
