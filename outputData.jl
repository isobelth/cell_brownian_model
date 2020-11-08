#
#  outputData.jl
#  cell-brownian-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module OutputData

using DelimitedFiles
using StaticArrays

@inline function outputData(pos::MMatrix,outfile::IOStream,t::Float64,tmax::Float64, age::MMatrix)
    outputdata = hcat(pos, age)
    writedlm(outfile, outputdata, ", " )
    #writedlm(outfile,pos,", ")
    flush(outfile)
    println("Simulating: $t/$tmax")

end

export outputData

end
