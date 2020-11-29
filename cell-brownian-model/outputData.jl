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

@inline function outputData(pos::MMatrix,outfile::IOStream,outfile2::IOStream, t::Float64,tmax::Float64, age::MMatrix, BM_pos::MMatrix)
    outputdata = hcat(pos, age)

    writedlm(outfile, outputdata, ", " )
    writedlm(outfile2, BM_pos, ",")
    flush(outfile)
    flush(outfile2)
    println("Simulating: $t/$tmax")

end

export outputData

end
