module JuliaBeagle

# Write your package code here.
using MCPhyloTree
using SharedArrays
using LinearAlgebra
using LoopVectorization

include("structs.jl")
include("algorithm.jl")
include("utils.jl")
end
