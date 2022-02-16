#Wrapper for variants of COMET/CRATER

module Multipoles

using COMET_None_mod
using COMET_CF_mod
using CRATER_None_mod
using CRATER_CF_mod
using CRATER_PCG_mod
using CRATER_PCG_CF_mod

using   Combinatorics,
        LinearAlgebra,
		DataFrames,
		Dates,
        LightGraphs

using CSV
using DataFrames
using Statistics

export	importData,
		COMET_None,
		COMET_CF,
		CRATER_None,
		CRATER_CF,
		CRATER_PCG,
		CRATER_PCG_CF

function importData(path)
	df = CSV.read(path, DataFrame,header = false)
    M = Matrix(df)
    println(size(M))
    C = cov(M)
    return C
end

function COMET_None(c::Array{F,2}, d::F, s::F, pair_lim::F) where {F<:AbstractFloat}
	M = COMET(c,d,s,pair_lim=pair_lim)
	return M
end

function COMET_CF(c::Array{F,2}, d::F, s::F, pair_lim::F) where {F<:AbstractFloat}
	M = COMET_CF_mod.COMET_CF(c,d,s,pair_lim=pair_lim)
	return M
end

function CRATER_None(c::Array{F,2}, d::F, s::F, pair_lim::F) where {F<:AbstractFloat}
	M = CRATER_None_mod.CRATER_None(c,d,s,pair_lim=pair_lim)
	return M
end

function CRATER_CF(c::Array{F,2}, d::F, s::F, pair_lim::F) where {F<:AbstractFloat}
	M = CRATER_CF_mod.CRATER_CF(c,d,s,pair_lim=pair_lim)
	return M
end

function CRATER_PCG(c::Array{F,2}, d::F, s::F, pair_lim::F) where {F<:AbstractFloat}
	M = CRATER_PCG_mod.CRATER_PCG(c,d,s,pair_lim=pair_lim)
	return M
end

function CRATER_PCG_CF(c::Array{F,2}, d::F, s::F, pair_lim::F) where {F<:AbstractFloat}
	M = CRATER_PCG_CF_mod.CRATER_PCG_CF(c,d,s,pair_lim=pair_lim)
	return M
end

end
