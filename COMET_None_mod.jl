#push!(LOAD_PATH, "./")

module COMET_None_mod

export  findCandidates,
        getMultipole,
        filterMultipole,
        removeRedundancy,
        COMET

include("remove_redun_adv.jl")

using   Combinatorics,
        LinearAlgebra,
		CSV,
		Tables
        # LightGraphs,
        # ImportMacros

using Common: BronKerbosch, convertNode, getNeighbors, BronKerboschL1
#using MasterRunner

# data structure of a clique
struct poles_t{F<:AbstractFloat,I<:Signed}
    dept::F # dependence
    gain::F	# gain
    mem::Vector{I} # indices of member timeseries
end

function getDependence(mat; myset=[])
    myset == [] && (myset = collect(1:size(mat,1)))
    return 1 - eigvals(mat[myset, myset])[1]
end

function getGain(mat; myset=[], dep=-1.0)
    myset == [] && (myset = collect(1:size(mat,1)))
    len_sub = length(myset)-1

    (dep < 0) && (dep = getDependence(mat, myset=myset))
    dep_sub = [getDependence(mat, myset=subset) for subset in combinations(myset, len_sub)]

    return dep - maximum(dep_sub)
end

function findCandidates(mat::Array{F, 2}; len_lim=3, pair_lim=1.0) where {F<:AbstractFloat}
	# number of timeseries
	num_var = size(mat, 1)

    mat_adj = zeros(Int, num_var*2, num_var*2)

    for i=1:num_var-1
        for j=i+1:num_var
            if abs(mat[i,j]) <= pair_lim
	            if mat[i,j] <= 0.0 # # add neg corr edges
	                mat_adj[i, j] = 1
	                mat_adj[j, i] = 1
	                mat_adj[i + num_var, j + num_var] = 1
	                mat_adj[j + num_var, i + num_var] = 1
	            else # add positive edges
	                mat_adj[i, j + num_var] = 1
	                mat_adj[j + num_var, i] = 1
	                mat_adj[i + num_var, j] = 1
	                mat_adj[j, i + num_var] = 1
				end
            end
		end
    end



	outfile = "C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\New\\ArchiveOutput\\COMET_AdjM.csv"
	#sopen(outfile, "w") do f
		#for i=1:num_var
			#for j=1:num_var
				#print(f,mat_adj[i,j])
			#end
			#print("\n")
		#end
		#print(f,mat_adj)
	#end
	CSV.write(outfile,  Tables.table(mat_adj), writeheader=false)

    cliques = Vector{Int64}[]
    # initialize Bron kerbosch on n nodes
    for node=1:num_var
        # println("Initialize with node $node")
        R = [node]
        X = collect(1:node-1)
        P = setdiff(getNeighbors(mat_adj, node), X)
		X = []
        BronKerbosch(mat_adj, R, P, X, cliques, min_size=len_lim)
		#BK L1 start here
		#mat_cor = abs.(mat)
		#mat_sq = mat .^ 2
		#L1_cur = 0.0
		#delta = 0.2
		#BronKerboschL1(mat_adj, mat_cor, R, P, X, cliques, L1_cur, delta, min_size=3, idx_lim=num_var, flag=true)
		#BK L1 end here
    end

	# G = []

    # numcliques = length(cliques)

    # canFlag = falses(numcliques)

    # for i=1:numcliques
    #     clique = cliques[i]
    #     minbefore = minimum(clique)
    #     clique[clique .> num_var] .-= num_var
    #     minafter = minimum(clique)
    #     minbefore == minafter ? canFlag[i] = true : canFlag[i] = false
    # end

    # deleteat!(cliques, findall(canFlag))

	# canFlag = []
	#println(findArrayInClique([78, 31, 89],cliques))
	return cliques
end

function getL1(x::Vector)
    return sum(abs.(x))
end

function getL1(mat::Array{Float64,2}, gr::Vector{Int})
    sz = size(mat,1)
    vec = [mat[row,col] for (row, col) in combinations(gr,2)]
    return getL1(vec)
end

function getL1Threshold(sz, delta)
    n_pair = Int(sz*(sz-1)/2)
    return n_pair*delta
end

function checkL1(mat, gr, delta)
    sz = length(gr)
    return getL1(mat, gr) >= getL1Threshold(sz, delta)
end

function getL2Sq(mat_sq, gr)
    sz = size(mat_sq,1)
    return sum(mat_sq[row,col] for (row,col) in combinations(gr,2))
end

function getL2boundSqGain(delta; sz=0)
    r = delta
    n_pair = sz*(sz-1)/2
    return n_pair*r^2
end

function getL2boundSqDependence(sigma; sz=0)
    r = sigma/(sz-1)
    n_pair = sz*(sz-1)/2
    return n_pair*r^2
end

function getL2boundSq(sz, delta, sigma)
    L2_delta = getL2boundSqGain(delta, sz=sz)
    L2_sigma = getL2boundSqDependence(sigma, sz=sz)
    return max(L2_delta, L2_sigma)
end

function checkL2(mat_sq, gr, sigma)
    return getL2Sq(mat_sq, gr) >= getL2boundSqDependence(sigma, sz=length(gr))
end

function getMultipole( mat, mat_sq, groups,
                        delta, sigma)
    subs = Vector{poles_t}()
    for gr in groups
        #pass_L2 = checkL2(mat_sq, gr, sigma)
        #!pass_L2 && (continue)

        dep = getDependence(mat, myset=gr)
        # skip if dep less than threshold
        dep < sigma && (continue)

        #pass_L1 = checkL1(mat, gr, delta)
        gain = getGain(mat, myset=gr, dep=dep)
        if gain >= delta # multipole: save and continue
            push!(subs, poles_t(dep, gain, gr))
            continue
        end


        # not a multipole, check all subsets that have the first 2 vertices
        sz = length(gr)
        sz == 3 && (continue)
        upperbound = Int(floor(1 + 1/delta))
        sublen = minimum([sz - 1, upperbound])

        while sublen >= 3
            for gr_sub in combinations(gr, sublen)
                #pass_L1 = checkL1(mat, gr_sub, delta)
                #pass_L2 = checkL2(mat_sq, gr_sub, sigma)

                #!(pass_L1 && pass_L2) && (continue)

                dep_sub = getDependence(mat, myset=gr_sub)
                dep_sub < sigma && (continue)
                gain_sub = getGain( mat, myset=gr_sub, dep=dep_sub)
                if gain_sub >= delta
                    push!(subs, poles_t(dep_sub, gain_sub, gr_sub))
                end
            end
            sublen -= 1
        end

    end

    return subs
end

# function filterMultipoles(M::Vector{poles_t})
#     # sort M
#     # sort!(M, by=x->length(x.mem), rev=true)

#     flag = falses(length(M))
# 	numpole = length(M)

#     for x=numpole:-1:2
#         flag[x] && (continue)
#         for i=x-1:-1:1
#             flag[i] && (continue)
#             if issubset(M[x].mem, M[i].mem)
#                 flag[x] = true
#                 break
#             elseif issubset(M[i].mem, M[x].mem)
#                 flag[i] = true
#             end
#         end
#     end

#     deleteat!(M, findall(flag))

#     return nothing
# end

function COMET(mat::Array{F,2}, delta::F, sigma::F; len_lim=3, pair_lim=1.0) where {F<:AbstractFloat}

    sz = size(mat, 1)

    candidates = findCandidates(mat, len_lim=len_lim, pair_lim=pair_lim)
	#lc = size(candidates)
	#println("Candidates: $lc")

	#outfile = "C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\New\\ArchiveOutput\\COMET_Candidates.txt"
	#open(outfile, "w") do f
	#	for i in eachindex(candidates)
	#		println(f,candidates[i])
	#	end
	#nd


    mat_sq = mat .^ 2

    M = getMultipole(mat,mat_sq, candidates, delta, sigma)
	#lm = size(M)
	#println("Multipoles: $lm")

	#outfile = "C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\New\\ArchiveOutput\\COMET_Multipoles.txt"
	#open(outfile, "w") do f
	#	for i in eachindex(M)
	#		println(f,M[i].mem)
	#	end
	#end
    # filterMultipoles(M)

    # return M
	M1 = M

	M1 = removeRedundancy(M1, sz)
	#M1 = removeDupes(M1)
	#lm1 = size(M1)
	#println("Max Multipoles: $lm1")
	return M1
	#outfile = "C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\New\\ArchiveOutput\\COMET_MultipolesMax.txt"
	#open(outfile, "w") do f
	#	for i in eachindex(M1)
	#		println(f,M1[i].mem)
	#	end
	#end
    #return M1
end

end
