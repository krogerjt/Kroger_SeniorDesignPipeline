# Similar to algorithm 22 except that L1 bound from adding two nodes
# is checked and used to break the loop as in algorithm 2


#push!(LOAD_PATH, "./")

module CRATER_PCG_CF_mod

export CRATER_PCG_CF

using Combinatorics, LinearAlgebra,CSV,Tables

# add common functions
using Common

# Wrapper for BronKerboschL1
# To call everytime after adding 1 edge to the graph
function findCliques(mat, mat_cor, delta;
                    has_edge=[], min_size=3,
                    output=[], idx_lim=0)
    R = []
    X = []
    P = collect(1:size(mat,1))
    idx_lim == 0 && (idx_lim = size(mat_cor, 1))
    # start fresh, no conditions on cliques to be found
    if has_edge==[]
        R = []
        X = []
        L1_cur = 0
    # cliques needs to have nodes in has_edge
    # initialize sets for recursion accordingly
    else
        R = has_edge[:]
        X = []
        sz = size(mat_cor,1)
        L1_cur = mat_cor[convertNode(sz, has_edge[1]), convertNode(sz, has_edge[2])]
        P = intersect(P, getNeighbors(mat, has_edge[1]), getNeighbors(mat, has_edge[2]))
    end

    BronKerboschL1(mat, mat_cor, R, P, X, output, L1_cur, delta, min_size=min_size, idx_lim=idx_lim, flag=true)

    return nothing
end

function CRATER_PCG_CF(mat, delta, sigma; pair_lim=1)
    sz = size(mat,1)
    # get elements and indices
    pairs = [[mat[row,col], (row, col)] for col=2:sz for row=1:col-1]
    # sort pairs on the abs. value
    sort!(pairs, by=x->abs(x[1]), rev=true)

    # index of starting pair
    if pair_lim == 0
        idx_start=1
    else
        idx_start = findfirst(x->abs(x[1])<=pair_lim, pairs)
    end

    # initialize empty adjacency matrix
    mat_adj = zeros(Int, 2*sz, 2*sz)

    # mat abs
    mat_cor = abs.(mat)
    mat_sq = mat .^ 2

    # initialize clique set
    cliques = []
    candidates_JK = []
    for idx=idx_start:length(pairs)
        pair = pairs[idx]
        # get val
        val = pair[1]

        # get indices
        row = pair[2][1]
        col = pair[2][2]
        has_edge = addEdge(mat_adj, val, row, col)

        # find clique containing the new edge
        cliques_new = []
        findCliques(mat_adj, mat_cor, delta, has_edge=has_edge, min_size=3, output=cliques_new, idx_lim=sz)

        cliques_new == [] && (continue)
        append!(candidates_JK, sort(cliques_new))
        # get multipoles from cliques
        cliques_new = assessMultipolesWithBounds(mat, mat_sq, cliques_new, delta=delta, sigma=sigma)

        cliques_new == [] && (continue)

        # # remove subsets
        # removeSubset!(cliques, cliques_new)

        # println(cliques_new)
        append!(cliques, cliques_new)
    end

	#outfile = "C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\New\\ArchiveOutput\\CRATER_AdjM.csv"
	#CSV.write(outfile,  Tables.table(mat_adj), writeheader=false)
    #outfile = "C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\New\\ArchiveOutput\\CRATER_Candidates.txt"
	#lc = size(candidates_JK)
	#println("Candidates: $lc")
	#lc = size(cliques)
	#println("Multipoles: $lc")
	#open(outfile, "w") do f
	#	for i in eachindex(candidates_JK)
	#		println(f,candidates_JK[i])
	#	end
	#end
	#outfile = "C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\New\\ArchiveOutput\\CRATER_Multipoles.txt"
	#open(outfile, "w") do f
	#	for i in eachindex(cliques)
	#		println(f,cliques[i].mem)
	#	end
	#end

	M1 = cliques
	#lc = size(M1)
	#println("Max Multipoles: $lc")
	M1 = removeRedundancy(M1, sz)
	#M1 = removeDupes(M1)
	#outfile = "C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\New\\ArchiveOutput\\CRATER_MultipolesMax.txt"
	#open(outfile, "w") do f
	#	for i in eachindex(M1)
	#		println(f,M1[i].mem)
	#	end
	#end

    return M1
end

end
