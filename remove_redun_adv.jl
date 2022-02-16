 #struct poles_t{F<:AbstractFloat,I<:Signed}
#     dept::F # dependence
#     gain::F	# gain
#     mem::Vector{I} # indices of member timeseries
 #end


mutable struct node_t
    id::Int # index of timeseries
    s::Int # support
    child::Vector{node_t} # child nodes
    dept::Float64
    gain::Float64
    function node_t()
        return new(0, 0, node_t[], -1.0, -1.0)
    end
    function node_t(x)
        return new(x, 0, node_t[], -1.0, -1.0)
    end
    function node_t(x,y,z)
        return new(x,y,z, -1.0, -1.0)
    end
end

function setQuality!(node; dept=-1.0, gain=-1.0)
    node.dept = dept
    node.gain = gain
    return nothing
end

function addChild!(node, childID)
    push!(node.child, node_t(childID))
    return nothing
end

function sortFreq(n, poles; rev=true)
    sup = zeros(Int, n)
    for x in poles
        series = x.mem
        sup[series] = sup[series] .+ 1
    end
    order = sortperm(sup,rev=rev)
    return [findfirst(order .== x) for x=1:n]
end

function sortLabels!(labels, order)
    sort!(labels, by=x->order[x])
    return nothing
end

# use ascending order of freq
function sortChildren!(node, order)
    if length(node.child) > 1
        sort!(node.child, by=x->order[x.id])
    end
    return nothing
end

function compareOrder(order, x, y, op)
    return op(order[x],order[y])
end

function updateNode!(node, series; idx_next=0, dept=-1.0, gain=-1.0)
    node.s = node.s + 1
    # failsafe
    if (idx_next ==  0)
        return nothing
    end
    # reaching leaf node
    if idx_next > length(series)
        setQuality!(node, dept=dept, gain=gain)
        return nothing
    end

    childID = series[idx_next]
    for child in node.child
        if child.id == childID
            updateNode!(child, series, idx_next=idx_next+1, dept=dept, gain=gain)
            return nothing
        end
    end

    # no child
    addChild!(node, childID)
    updateNode!(node.child[end], series, idx_next=idx_next+1, dept=dept, gain=gain)
    return nothing
end

# use descending order of freq
function addPole2Tree!(tree, order_rev, series; dept=-1.0, gain=-1.0)
    # sort series
    sortLabels!(series, order_rev)
    updateNode!(tree[series[1]], series, idx_next=2, dept=dept, gain=gain)
    return nothing
end

# use descending order of freq
function buildTree(n_region, poles, order_rev)
    tree = [node_t(x) for x=1:n_region]
    for pole in poles
        series = pole.mem
        # update tree with new pole
        addPole2Tree!(tree, order_rev, series, dept=pole.dept, gain=pole.gain)
    end
    return tree
end

# use ascending order of freq
function getPath2Leaf(node, path, order)
    push!(path, node.id)
    if node.child == []
        if length(path) >= 3
            return [poles_t(node.dept, node.gain, Int.(path[:]))]
        else
            return []
        end
    end
    output = []
    sortChildren!(node, order)
    for child in node.child
        append!(output, getPath2Leaf(child, path[:], order))
    end

    return output
end

function getSets(tree, order)
    output = []
    tree_order = sortperm(order)
    for node in tree[tree_order]
        node.child == [] && (continue)
        append!(output, getPath2Leaf(node, [], order))
    end

    return output
end

mutable struct nodeLex_t
    b::Int
    e::Int
    j::Int
    d::Int
    m::Int
    child_left::Union{nodeLex_t, Nothing}
    child_right::Union{nodeLex_t, Nothing}

    function nodeLex_t()
        return new(0,0,0,0,0,nothing, nothing)
    end

    function nodeLex_t(b,e,j,d,m)
        return new(b,e,j,d,m, nothing, nothing)
    end
end

function getNextItem(S::Vector{Int}, j::Int, item_id::Int, order)
    S_ordered = [compareOrder(order, x, item_id, >=) for x in S[j+1:end]]
    idx_find = findfirst(S_ordered)
    isnothing(idx_find) && (return nothing)
    return idx_find + j
end

# find the index of the last group within b:e such that the d_th item equals item_id
function getNextEndRange(D::Vector{Vector{Int}}, b::Int, e::Int, item_id::Int, d::Int, order)
    # get ids of items at position d
    items_d = [x[d] for x in D[b:e]]
    # convert into order
    items_ordered = [compareOrder(order, x, item_id, ==) for x in items_d]
    # find the last set in D[b:e] where item at d equals item_id
    idx_next = findlast(items_ordered)
    isnothing(idx_next) && (return nothing)
    return  idx_next + b - 1
end

# find the first sets in D[b:e] whose item at position d is different from item_id
function getNextBeginRange(D::Vector{Vector{Int}}, b::Int, e::Int, item_id::Int, d::Int, order)
    # get ids of items at position d
    items_d = [x[d] for x in D[b:e]]
    # compare with item_id
    items_ordered = [compareOrder(order, x, item_id, >=) for x in items_d]

    # find the last set in D[b:e] where item at d >= item_id
    idx_next = findfirst(items_ordered)
    isnothing(idx_next) && (return e+1)
    return idx_next + b - 1
end

# check for subsets of S in D[b:e] starting from position j
# all sets in D[b:e] need to share a prefix of length d with D[b]
function markSubsumed(D, b, e, S, j, d, parent, flags, order)
    # println(tab*"b=", b, "\te=", e, "\tj=", j, "\td=", d)
    b > e && (return nothing)
    b = findnext(flags, b)
    if isnothing(b) || (b > e)
        # println(tab*"No b")
        return nothing
    end

    len_S = length(S)

    # S[j] < D[b][d+1]
    if compareOrder(order, S[j], D[b][d+1], <)
        # println("\tLess")
        j = getNextItem(S, j, D[b][d+1], order)
        if isnothing(j) # nothing to do for this node
            # assign a big value to m so that the node will not be reused
            parent.m = 999999
            return nothing
        end
        # println(tab*"New j = ", j)
    end
    parent.m = j

    # S[j] = D[b][d+1]
    if compareOrder(order, S[j], D[b][d+1], ==) # match

        e_new = getNextEndRange(D, b, e, S[j], d+1, order)
        # println(tab*"New e = ", e_new)
        # find and mark all obvious subsets
        if len_S > d + 1
            # println("\tD[b]=", length(D[b]), "\td+1=", d+1)
            while (b <= e_new) && (length(D[b]) == d + 1)
                # println("\tD[b] = d + 1 = ", d+1)
                flags[b] = false
                # println("\tFlag[b] = ", flags[b])
                # println(tab*"\tFound : ", b)
                b += 1
            end
        end

        if (j + 1 <= len_S) && (b <= e_new)
            # left part
            parent.child_left = nodeLex_t(b, e_new, j+1, d+1, j+1)
            # println(tab*"Go left")
            markSubsumed(D, b, e_new, S, j+1, d+1, parent.child_left, flags, order)
        end

        # right part
        b = e_new + 1
        parent.child_right = nodeLex_t(b, e, j, d, j)
        # println(tab*"Go right")
        markSubsumed(D, b, e, S, j, d, parent.child_right, flags, order)

    else
        # S[j] > D[b][d+1]
        # chop first few sets until item at index d is different
        b = getNextBeginRange(D, b, e, S[j], d+1, order)
        # println(tab*"New b: ", b)
        # println(tab*"Pruning head")
        parent.child_right = nodeLex_t(b, e, j, d, j)
        markSubsumed(D, b, e, S, j, d, parent.child_right, flags, order)
    end

    return nothing
end

function checkSet(D::Vector{Vector{Int}}, i::Int, p::Int, v::nodeLex_t, order, flags)

        # # modify the node by updating b, removing children and start new
        # v.child_left = nothing
        # v.child_right = nothing
        # markSubsumed(D, v.b, v.e, D[i], v.j, v.d, v, flags, order)

    if v.m > p
        b = max(v.b, i+1)
        # modify the node by updating b, removing children and start new
        v.b = b
        v.child_left = nothing
        v.child_right = nothing
        # println(tab*"Do new for m = $(v.m)\tp = $p")
        markSubsumed(D, v.b, v.e, D[i], v.j, v.d, v, flags, order)
        return nothing
    end



    if !isnothing(v.child_left)
        checkSet(D, i, p, v.child_left, order, flags)
    end

    if !isnothing(v.child_right)
        checkSet(D, i, p, v.child_right, order, flags)
    end
    return nothing
end

function removeNonPrefixedSubsets(D, order)
    n_set = length(D)
    flags = trues(n_set)
    v = nothing
    S = nothing
    for i=1:n_set-1
        flags[i] || (continue)
        if isnothing(v) # and also S
            b, e, j, d, m = 1, n_set, 1, 0, 1
            v = nodeLex_t(b, e, j, d, m)
            p = 0
        else
            # find p
            p = 0
            for pos=1:min(length(S), length(D[i]))
                if compareOrder(order, S[pos], D[i][pos], ==)
                    p += 1
                else
                    break
                end
            end
        end

        checkSet(D, i, p, v, order, flags)
        S = D[i]
        # break
    end

    return flags
end

# get series from all multipoles
function getSeries(poles)
    return [pole.mem for pole in poles]
end

function findUniqueMaximal(M, n_region)
    # get timeseries indices in descending order of support
    order_rev = sortFreq(n_region, M, rev=true)
    # println("Got order_rev")
    # println("Got order")
    # build the FP Growth tree
    tree = buildTree(n_region, M, order_rev)
    # println("Built tree")
    # get multipoles
    poles = getSets(tree, order_rev)
    # println("Got sets from FP tree")
    # get only series
    poles_series = getSeries(poles)
    # println("\tN_pole before last step: ", length(poles_series))
    # println("\nN_nofilter: ", length(M), "\tN_semifilter: ", length(poles_series))
    flags = removeNonPrefixedSubsets(poles_series, order_rev)
    outPoles = []
    for i=1:length(poles_series)
        !(flags[i]) && (continue)
        push!(outPoles,poles[i])
    end

    return outPoles
end

function removeRedundancy(M, n_region)
    #ids = findall(findUniqueMaximal(M, n_region))
    poles = findUniqueMaximal(M, n_region)
    #println(ids)
    #println(id2)
    return poles
end
