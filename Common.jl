module Common

export  poles_t,
        getDependence,
        getGain,
        getL1,
        getL1Threshold,
        checkL1,
        checkL2,
        getNeighbors,
        convertNode,
        toContinue,
        toStop,
        toRepeat,
        toTerminate,
        removeSubset!,
        addEdge,
        sortNode,
        BronKerbosch,
        BronKerboschL1,
        BronKerboschPivot,
        generateGraphMatrix,
        assessMultipoles,
        assessMultipolesWithBounds,
        assessMultipolesWithL1Bounds,
        removeRedundancy

using Combinatorics, LinearAlgebra
include("remove_redun_adv.jl")

# data structure of a clique
struct poles_t{F<:AbstractFloat,I<:Signed}
    dept::F # dependence
    gain::F	# gain
    mem::Vector{I} # indices of member timeseries
end

# convert node indices from scale 1:2m to scale 1:m
function convertNode(sz, node)
    node > sz ? (return node - sz) : (return node)
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

function getL1(x::Vector)
    return sum(abs.(x))
end

function getL1(mat::Array{Float64,2}, gr::Vector{Int})
    sz = size(mat,1)
    gr_conv = convertNode.(sz, gr)
    vec = [abs(mat[row,col]) for (row, col) in combinations(gr_conv,2)]
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
    gr_conv = convertNode.(sz, gr)
    return sum(mat_sq[row,col] for (row,col) in combinations(gr_conv,2))
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

function generateAdjacencyMat(mat, range_val)
    n_row, n_col = size(mat)
    val_min = range_val[1]
    val_max = range_val[2]
    output = [val_min < mat[row,col] < val_max ? 1 : 0 for col=1:n_col,row=1:n_row]
    output = output + output'
    return output
end

function getNeighbors(mat, node)
    return findall(mat[node, :] .== 1)
end

function getNeighborsBad(mat, node)
    N = findall(mat[node, :] .== 1)
    sz = Int((size(mat,1))/2)
    out = zeros(size(N[i]))
    for i in length(N)
        if N[i] > sz
            out[i] = Int(N[i] -sz)
        else
            out[i] = Int(N[i])
        end
    end
    return out
end

# update the sets to continue the recursion in finding cliques
function toContinue(mat, node, P, X, R)
    neighbors = getNeighbors(mat, node)
    R_temp = vcat(R, node)
    P_temp = intersect(P, neighbors)
    X_temp = intersect(X, neighbors)
    return P_temp, X_temp, R_temp
end

# update the sets to meet the stopping condition of the recursion in finding cliques
function toStop(X, R)
    return [], X, R
end

# maintain the sets to restart the recursion in finding cliques
function toRepeat(P, X, node_cur)
    # P_temp = P
    # X_temp = X
    # R_temp = node_cur
    return P[:], X[:], node_cur[:]
end

function toTerminate()
    return [], [1], []
end

# remove elements that are subsets of others
function removeSubset!(cliques)
    cliques == [] && (return nothing)
    n_clique = length(cliques)

    flags = falses(n_clique)

    for idx_left=1:n_clique
        flags[idx_left] && (continue)
        left = cliques[idx_left].mem
        for idx_right=1:n_clique
            idx_left == idx_right && (continue)
            flags[idx_right] && (continue)
            right = cliques[idx_right].mem
            if issubset(left, right)
                flags[idx_left]= true
                break
            end
        end
    end
    deleteat!(cliques, findall(flags))
    return nothing
end

# remove elements in cliques that are subsets of those in cliques_new
function removeSubset!(cliques, cliques_new)
    (cliques == [] || cliques_new == []) && (return nothing)
    n_clique = length(cliques)

    flags = falses(n_clique)

    for idx_old=1:n_clique
        old = cliques[idx_old].mem
        for new in cliques_new
            if issubset(old, new.mem)
                flags[idx_old]= true
                break
            end
        end
    end
    deleteat!(cliques, findall(flags))
    return nothing
end

# add one edge to the adjacency matrix in scale 1:2m according to the sign of the edge weight
function addEdge(mat::Array{Int,2}, val::Float64, row::Int, col::Int)
    sz = Int(size(mat,1)/2)
    if val <= 0.0
        # add to 1st half
        mat[row, col] = 1
        mat[col, row] = 1
        # add to 2nd half (mirror)
        mat[row+sz, col+sz] = 1
        mat[col+sz, row+sz] = 1
        has_edge = [row, col]
    else
        # add between 2 halves
        mat[row, col+sz] = 1
        mat[col+sz, row] = 1
        mat[col, row+sz] = 1
        mat[row+sz, col] = 1
        has_edge=[row, col+sz]
    end

    return has_edge
end

# generate empty adjacency matrix of scale 2m
function generateGraphMatrix(mat_cor)
    sz = size(mat_cor,1)
    mat = zeros(Int, 2*sz, 2*sz)
    for col=2:sz, row=1:col-1
        val = mat_cor[row,col]
        addEdge(mat, val, row, col)
    end

    return mat
end

# vanilla algorithm in finding cliques
function BronKerbosch(mat, R, P, X, output; idx_lim=0, min_size=3)
    if P == [] # no more common neighbor can be added to R
        if X == [] # is a maximum clique
            if length(R) >= min_size # meets the size threshold
                push!(output, [convertNode(idx_lim,node) for node in R])
            end
        end
        return true
    end

    idx_lim == 0 && (idx_lim = Int(size(mat,1)/2))
    # println(idx_lim)
    for p in P
        neighbors = getNeighbors(mat, p)
        R_temp = vcat(R, p)

        # println(R_temp[end])

        P_temp = intersect(P, neighbors)
        X_temp = intersect(X, neighbors)
        BronKerbosch(mat, R_temp, P_temp, X_temp, output, idx_lim=idx_lim, min_size=min_size)
        P = setdiff(P, p)
        push!(X, p)
    end

    return true
end

# compute the sum weights of edges between each node in nodes_new
# to all nodes in nodes_cur then sort nodes_new in the descending order
# Return the vec of sum weights for each node in nodes_new
function sortNode(mat, R, P)
    sz = size(mat, 1)
    # convert node idx back to scale 1:sz from scale 1:2*sz
    R_org = convertNode.(sz, R)
    P_org = convertNode.(sz, P)
    sum_cor = sum(mat[P_org, R_org], dims=2)[:]
    ids_sort = sortperm(sum_cor, rev=true)
    P[:] = P[ids_sort]

    return sum_cor[ids_sort]
end

function BronKerboschL1(mat, mat_cor,
                        R, P, X,
                        output, L1, delta;
                        idx_lim=0, min_size=3,
                        r_max=0.4, flag=true)

    if P == [] # no more common neighbor can be added to R
        if X == [] # is a maximum clique
            if length(R) >= min_size # meets the size threshold
                push!(output, [convertNode(idx_lim,node) for node in R])
                return true
            end
        end
        return false
    end

    # make a copy of X
    X_org = X[:]

    # get size of the current clique
    sz_cur = length(R)

    # flag that indicates whether a maximal clique is found from any child
    got_sth = false

    # compute the sum weights and sort P
    w_P = sortNode(mat_cor, R, P)

    # enumerate nodes in P and make recursive call
    to_break = false
    for (p, w_p) in zip(P, w_P)
        # generate sets to continue the recursion by default
        P_temp, X_temp, R_temp = toContinue(mat, p, P, X, R)

        # check L1 vs threshold
        # L1_next = getL1(mat_cor, vcat(R, p))
        L1_next = L1 + w_p

        if L1_next >= getL1Threshold(sz_cur+1, delta)
            flag_next = true
        else
            flag_next = false
            if !flag # two consecutive failures so terminate the recursion
                P_temp, X_temp, R_temp = toTerminate()
            else
                # check whether L1 bound on adding two nodes also passes
                L1_bound = L1_next + w_p + r_max
                if L1_bound < getL1Threshold(sz_cur+2, delta)
                    P_temp, X_temp, R_temp = toStop(X, R)
                    to_break = true
                end
            end
        end

        got_sth = got_sth | BronKerboschL1( mat, mat_cor,
                                R_temp, P_temp, X_temp,
                                output, L1_next, delta,
                                idx_lim=idx_lim, min_size=min_size, flag=flag_next)
        to_break && (break)
        # update P and X before enumerating the next node in P
        if flag_next
            P = setdiff(P, p)
            push!(X, p)
        end
    end

    # check whether no maximal clique was found in the for loop
    if !got_sth
        if flag # R satisfies L1 threshold so output R as a maximal clique
            P_temp, X_temp, R_temp = toStop(X_org, R)
        else # otherwise output nothing and terminate
            P_temp, X_temp, R_temp = toTerminate()
        end
        return BronKerboschL1(  mat, mat_cor,
                    R_temp, P_temp, X_temp,
                    output, L1, delta,
                    idx_lim=idx_lim, min_size=min_size, flag=flag)
    else
        return true
    end
end

# implement the Bron Kerbosch variant by choosing a pivot among P
# and only recursively call on the pivot and its non-neighbors within P
function BronKerboschPivot(mat, R, P, X, output; idx_lim=0, min_size=3)
    if (P == []) && (X == [])
        if length(R) >= min_size
            push!(output, [convertNode(idx_lim,node) for node in R])
            # println(R)
        end
        return nothing
    end
    idx_lim == 0 && (idx_lim = size(mat,1))
    # println(idx_lim)

    # randomly choose a pivot
    pivot = rand(vcat(P, X))
    # iterate all nodes in P except neighbors to pivot
    for node in setdiff(P, getNeighbors(mat, pivot))
        neighbors = getNeighbors(mat, node)
        R_temp = vcat(R, node)

        # println(R_temp[end])

        P_temp = intersect(P, neighbors)
        X_temp = intersect(X, neighbors)
        BronKerbosch(mat, R_temp, P_temp, X_temp, output, idx_lim=idx_lim, min_size=min_size)
        P = setdiff(P, [node])
        push!(X, node)
    end

    return nothing
end

# improve BronKerboschPivot2 by adding a set for nodes_pivot
# only enumerate from nodes_pivot instead of P
# Still have flaws to be fixed
function BronKerboschPivot2(mat, R, P, X, nodes_pivot, output; idx_lim=0, min_size=3)
    if (P == []) && (X == [])
        if length(R) >= min_size
            push!(output, [convertNode(idx_lim,node) for node in R])
            # println(R)
        end
        return nothing
    end
    idx_lim == 0 && (idx_lim = size(mat,1))
    # println(idx_lim)

    # randomly choose a pivot
    nodes_pivot == [] && (nodes_pivot = P[:])
    pivot = rand(nodes_pivot)


    neighbors = getNeighbors(mat, pivot)

    # add pivot to R and call function on neighbors of pivot
    R_temp = vcat(R, pivot)
    P_temp = intersect(P, neighbors)
    X_temp = intersect(X, neighbors)
    nodes_pivot_temp = []
    BronKerboschPivot2(mat, R_temp[:], P_temp[:], X_temp[:], nodes_pivot_temp[:], output, idx_lim=idx_lim, min_size=min_size)

    # recall the function on the non-neighbors of pivot
    R_temp = R
    # nothing is added to R, so only need to remove pivot from P
    P = setdiff(P, pivot)
    cluster_pivot = vcat(pivot, neighbors)
    push!(X, pivot)
    nodes_pivot = setdiff(nodes_pivot, cluster_pivot)
    BronKerboschPivot2(mat, R_temp[:], P[:], X[:], nodes_pivot[:], output, idx_lim=idx_lim, min_size=min_size)

    return nothing
end

function assessMultipoles(  mat, groups; delta=0.1, sigma=0.3)
    subs = []
    for gr in groups
        dep = getDependence(mat, myset=gr)
        # skip if dep less than threshold
        dep < sigma && (continue)

        gain = getGain(mat, myset=gr, dep=dep)
        if gain >= delta # multipole: save and continue
            push!(subs, poles_t(dep, gain, gr))
            continue
        end

        # not a multipole, check all subsets that have the first 2 vertices
        sz = length(gr)
        sz == 3 && (continue)
        pair = gr[1:2]
        for sz_sub=sz-3:-1:1
            for others in combinations(gr[3:end], sz_sub)
                gr_sub = vcat(pair, others)
                dep_sub = getDependence(mat,
                                myset=gr_sub)
                dep_sub < sigma && (continue)
                gain_sub = getGain( mat,
                            myset=gr_sub,
                            dep=dep_sub)
                if gain_sub >= delta
                    push!(subs, poles_t(dep_sub,
                                gain_sub,
                                gr_sub))
                end
            end
        end
    end

    # remove new multipoles found as subsets of cliques
    removeSubset!(subs)
    return subs
end

function assessMultipolesWithBounds(  mat, mat_sq, groups;
                            delta=0.1, sigma=0.3)
    # n_gr = length(groups)
    # flags = trues(n_gr)
    subs = []
    for gr in groups
        pass_L2 = checkL2(mat_sq, gr, sigma)
        pass_L2 || (continue)

        dep = getDependence(mat, myset=gr)
        # skip if dep less than threshold
        dep < sigma && (continue)

        pass_L1 = checkL1(mat, gr, delta)
        if pass_L1
            gain = getGain(mat, myset=gr, dep=dep)
            if gain >= delta # multipole: save and continue
                push!(subs, poles_t(dep, gain, gr))
                continue
            end
        end

        # not a multipole, check all subsets that have the first 2 vertices
        sz = length(gr)
        sz == 3 && (continue)
        pair = gr[1:2]
        for sz_sub=sz-3:-1:1
            for others in combinations(gr[3:end], sz_sub)
                gr_sub = vcat(pair, others)

                pass_L1 = checkL1(mat, gr_sub, delta)
                pass_L2 = checkL2(mat_sq, gr_sub, sigma)

                !(pass_L1 && pass_L2) && (continue)

                dep_sub = getDependence(mat,
                                        myset=gr_sub)
                dep_sub < sigma && (continue)
                gain_sub = getGain( mat,
                                    myset=gr_sub,
                                    dep=dep_sub)
                if gain_sub >= delta
                    push!(subs, poles_t(dep_sub,
                                    gain_sub,
                                    gr_sub))
                end
            end
        end
    end

    # remove new multipoles found as subsets of cliques
    removeSubset!(subs)
    return subs
end

function assessMultipolesWithL1Bounds(  mat, groups;
                                        delta=0.1, sigma=0.3)
    # n_gr = length(groups)
    # flags = trues(n_gr)
    subs = []
    for gr in groups
        dep = getDependence(mat, myset=gr)
        # skip if dep less than threshold
        dep < sigma && (continue)

        pass_L1 = checkL1(mat, gr, delta)
        if pass_L1
            gain = getGain(mat, myset=gr, dep=dep)
            if gain >= delta # multipole: save and continue
                push!(subs, poles_t(dep, gain, gr))
                continue
            end
        end

        # not a multipole, check all subsets that have the first 2 vertices
        sz = length(gr)
        sz == 3 && (continue)
        pair = gr[1:2]
        for sz_sub=sz-3:-1:1
            for others in combinations(gr[3:end], sz_sub)
                gr_sub = vcat(pair, others)

                pass_L1 = checkL1(mat, gr_sub, delta)

                !pass_L1 && (continue)

                dep_sub = getDependence(mat,
                                myset=gr_sub)
                dep_sub < sigma && (continue)
                gain_sub = getGain( mat,
                            myset=gr_sub,
                            dep=dep_sub)
                if gain_sub >= delta
                    push!(subs, poles_t(dep_sub,
                                gain_sub,
                                gr_sub))
                end
            end
        end
    end

    # remove new multipoles found as subsets of cliques
    removeSubset!(subs)
    return subs
end

end
