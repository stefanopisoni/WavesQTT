module WavesQTT

using ITensors
using ITensorMPS
using ITensorQTT
using LinearAlgebra #Diagonal
using FFTW

include("splitstep.jl")

export qtt_splitsteps

end
