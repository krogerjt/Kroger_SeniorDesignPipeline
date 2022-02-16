module MasterRunner

using COMET_None_mod
using COMET_CF_mod
using CRATER_None_mod
using CRATER_CF_mod
using CRATER_PCG_mod
using CRATER_PCG_CF_mod

include("importCSV.jl")
include("remove_redun_adv.jl")

export 	Runner,importer, looper, Archive,poleToString,comparePoles,getGainStats,getDeptStats,sizeStats, Runner_mults, Lex,
		findPolesWithID,equals_PoleArray, findArrayInPoles,findArrayInClique,removeDupes,isInOrIsSubsetIn,findPolesWithID_DG, printAdjOfPole,adjPrintAll,importAdj,
		Pipeline, Pipeline_Mults

using   Combinatorics,
        LinearAlgebra,
		DataFrames,
		Dates,
        LightGraphs



function Runner(d::F, s::F,pair_lim::F) where {F<:AbstractFloat}
	carrier = DataFrame(type = Any[],time = Any[], multipoles = Any[])
	println("Starting runner\nStarting CSV Import")
	c = importCSV()
	time0 = time()
	println("Regular COMET")
	M = COMET(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("COMET done. Time = $time1 Mults = $MS")
	push!(carrier,("COMET Lex",time1,MS))

	time0 = time()
	println("COMET L1")
	#function COMET_CF(mat::Array{F,2}, delta::F, sigma::F; len_lim=3, pair_lim=1.0) where {F<:AbstractFloat}
	M = COMET_CF(c,d,s,len_lim=3.0,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("COMET L1 done. Time = $time1 Mults = $MS")
	push!(carrier,("COMET L1",time1,MS))

	time0 = time()
	println("CRATER Plain")
	#function CRATERPlain(mat, delta, sigma; pair_lim=1)
	M = CRATER_None(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER Plain done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER Plain",time1,MS))

	time0 = time()
	println("CRATER CF Only")
	#function CRATERCF(mat, delta, sigma; pair_lim=1)
	M = CRATER_CF(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER CF Only done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER CF Only",time1,MS))

	time0 = time()
	println("CRATER PCG Only")
	#function CRATERPCG(mat, delta, sigma; pair_lim=1)
	M = CRATER_PCG(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER PCG Only done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER PCG Only",time1,MS))

	time0 = time()
	println("CRATER Full")
	#function CRATER(mat, delta, sigma; pair_lim=1)
	M = CRATER_PCG_CF(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER Full done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER Full",time1,MS))

	#delete!(carrier,["filler",0.0,0])
	#CSV.write("C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\outputMasterRunner1116 $d $s $pair_lim.csv",c2,header = true)

	return carrier
end

function Runner_mults(d::F, s::F,pair_lim::F) where {F<:AbstractFloat}
	carrier = DataFrame(type = Any[],time = Any[], multipoles = Any[])
	println("Starting runner\nStarting CSV Import")
	c = importCSV()
	time0 = time()
	println("Regular COMET")
	M = COMET(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("COMET done. Time = $time1 Mults = $MS")
	push!(carrier,("COMET Lex",time1,MS))
	m1 = M

	time0 = time()
	println("COMET L1")
	#function COMET_CF(mat::Array{F,2}, delta::F, sigma::F; len_lim=3, pair_lim=1.0) where {F<:AbstractFloat}
	M = COMET_CF(c,d,s,len_lim=3.0,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("COMET L1 done. Time = $time1 Mults = $MS")
	push!(carrier,("COMET L1",time1,MS))
	m2 = M

	time0 = time()
	println("CRATER Plain")
	#function CRATERPlain(mat, delta, sigma; pair_lim=1)
	M = CRATER_None(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER Plain done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER Plain",time1,MS))
	m3 = M

	time0 = time()
	println("CRATER CF Only")
	#function CRATERCF(mat, delta, sigma; pair_lim=1)
	M = CRATER_CF(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER CF Only done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER CF Only",time1,MS))
	m4 = M

	time0 = time()
	println("CRATER PCG Only")
	#function CRATERPCG(mat, delta, sigma; pair_lim=1)
	M = CRATER_PCG(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER PCG Only done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER PCG Only",time1,MS))
	m5 = M

	time0 = time()
	println("CRATER Full")
	#function CRATER(mat, delta, sigma; pair_lim=1)
	M = CRATER_PCG_CF(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER Full done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER Full",time1,MS))
	m6 = M

	CSV.write("C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\outputMasterRunner1116 $d $s $pair_lim.csv",c2,header = true)

	return carrier,m1,m2,m3,m4,m5,m6
end

function Pipeline(c::Array{F,2}, d::F, s::F,pair_lim::F) where {F<:AbstractFloat}
	#Run Pipeline with custom cov matrix and parameters
	carrier = DataFrame(type = Any[],time = Any[], multipoles = Any[])
	time0 = time()
	println("Regular COMET")
	M = COMET(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("COMET done. Time = $time1 Mults = $MS")
	push!(carrier,("COMET Lex",time1,MS))

	time0 = time()
	println("COMET L1")
	#function COMET_CF(mat::Array{F,2}, delta::F, sigma::F; len_lim=3, pair_lim=1.0) where {F<:AbstractFloat}
	M = COMET_CF(c,d,s,len_lim=3.0,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("COMET L1 done. Time = $time1 Mults = $MS")
	push!(carrier,("COMET L1",time1,MS))

	time0 = time()
	println("CRATER Plain")
	#function CRATERPlain(mat, delta, sigma; pair_lim=1)
	M = CRATER_None(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER Plain done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER Plain",time1,MS))

	time0 = time()
	println("CRATER CF Only")
	#function CRATERCF(mat, delta, sigma; pair_lim=1)
	M = CRATER_CF(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER CF Only done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER CF Only",time1,MS))

	time0 = time()
	println("CRATER PCG Only")
	#function CRATERPCG(mat, delta, sigma; pair_lim=1)
	M = CRATER_PCG(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER PCG Only done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER PCG Only",time1,MS))

	time0 = time()
	println("CRATER Full")
	#function CRATER(mat, delta, sigma; pair_lim=1)
	M = CRATER_PCG_CF(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER Full done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER Full",time1,MS))

	#CSV.write("C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\outputMasterRunner1116 $d $s $pair_lim.csv",c2,header = true)

	return carrier
end

function Pipeline_Mults(c::Array{F,2}, d::F, s::F,pair_lim::F) where {F<:AbstractFloat}
	carrier = DataFrame(type = Any[],time = Any[], multipoles = Any[])
	time0 = time()
	println("Regular COMET")
	M = COMET(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("COMET done. Time = $time1 Mults = $MS")
	push!(carrier,("COMET Lex",time1,MS))
	m1 = M

	time0 = time()
	println("COMET L1")
	#function COMET_CF(mat::Array{F,2}, delta::F, sigma::F; len_lim=3, pair_lim=1.0) where {F<:AbstractFloat}
	M = COMET_CF(c,d,s,len_lim=3.0,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("COMET L1 done. Time = $time1 Mults = $MS")
	push!(carrier,("COMET L1",time1,MS))
	m2 = M

	time0 = time()
	println("CRATER Plain")
	#function CRATERPlain(mat, delta, sigma; pair_lim=1)
	M = CRATER_None(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER Plain done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER Plain",time1,MS))
	m3 = M

	time0 = time()
	println("CRATER CF Only")
	#function CRATERCF(mat, delta, sigma; pair_lim=1)
	M = CRATER_CF(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER CF Only done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER CF Only",time1,MS))
	m4 = M

	time0 = time()
	println("CRATER PCG Only")
	#function CRATERPCG(mat, delta, sigma; pair_lim=1)
	M = CRATER_PCG(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER PCG Only done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER PCG Only",time1,MS))
	m5 = M

	time0 = time()
	println("CRATER Full")
	#function CRATER(mat, delta, sigma; pair_lim=1)
	M = CRATER_PCG_CF(c,d,s,pair_lim=pair_lim)
	MS = length(M)
	time1 = time() - time0
	println("CRATER Full done. Time = $time1 Mults = $MS")
	push!(carrier,("CRATER Full",time1,MS))
	m6 = M

	#CSV.write("C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\outputMasterRunner1116 $d $s $pair_lim.csv",c2,header = true)

	return carrier,m1,m2,m3,m4,m5,m6
end

function importer()
	carrier = DataFrame(type = ["filler"],time = [0.0], multipoles = [0])
	println("Starting runner\nStarting CSV Import")
	c = importCSV()
	return c
end

function looper()
	for miu_L in [1,0.95,0.9]
		for d_L in [0.2,0.21,0.19]
			for s_L in [0.5,0.51,0.49]
					d = d_L #Delta
					s = s_L #Sigma
					pair_lim = miu_L #miu
					println("Using parameters\n\tDelta = ",d,"\n\tSigma = ",s,"\n\tPL = ",pair_lim)
					Runner(d,s,pair_lim)
				end
			end
		end
end

function Archive()
	carrier = DataFrame(type = ["filler"],time = [0.0], multipoles = [0])
	println("Starting Archive runner\nStarting CSV Import")
	c = importCSV()
	d=0.2
	s = 0.6
	pair_lim = 1
	time0 = time()
	println("COMET L1")
	#function COMET_CF(mat::Array{F,2}, delta::F, sigma::F; len_lim=3, pair_lim=1.0) where {F<:AbstractFloat}
	M1 = COMET_CF(c,d,s,len_lim=0.0,pair_lim=pair_lim)
	if !isnothing(M1)
		MS1 = length(M1)
	else
		MS1 = 0
	end
	time1 = time() - time0
	println("COMET L1 done. Time = $time1 Mults = $MS1")
	push!(carrier,("COMET L1",time1,MS1))

	time0 = time()
	println("CRATER Full")
	#function CRATER(mat, delta, sigma; pair_lim=1)
	M2 = CRATER_PCG_CF(c,d,s,pair_lim=pair_lim)
	if !isnothing(M2)
		MS2 = length(M2)
	else
		MS2 = 0
	end
	time1 = time() - time0
	println("CRATER Full done. Time = $time1 Mults = $MS2")
	push!(carrier,("CRATER Full",time1,MS2))

	CSV.write("C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\outputMasterArchiveRunner1116 $d $s $pair_lim.csv",carrier,header = true)
	return M1, M2



end

function poleToString(M)

	ostring = []
	for poleID=1:length(M)
		arr = M[poleID].mem
		stringNew = ""
		for mem=1:length(arr)
			stringNew = stringNew * string(arr[mem])*","
		end
		ostring = stringNew
	end
	return ostring
end

function comparePoles(M1,M2)
	flag = true
	outID = []							#initalize the flag boolean
	for smallpoleID=1:length(M1)		#iterate through each member of the multipoles with less members
		flag = false					#set the flag to false before checking against the larger set
		sP = sort(M1[smallpoleID].mem)	#sort the current searched member of the smaller pole in number order. ie [3,2,1] becomes [1,2,3]
		for bigpoleID=1:length(M2)		#iterate through each member of the larger set for each member of the smaller set
			bP = sort(M2[bigpoleID].mem)#same sorting logic as the smaller poles for the current iteratation of the larger pole
			if issubset(sP,bP)				#straight comparison of the two sorted arrays
				flag = true				#if the two arrays are equal, the smaller pole is a member of the larger poles
			end
		end
		if flag == false									#if a small pole array made it through the iteration of the larger poles without triggering the flag bool to true, the smaller collection
			outID = append!(outID,smallpoleID)				#	must not be a subset of the larger group, therefore return false.
		end												#	could be modified to return a list of non-subset poles and their dep/gain
	end
	return findPolesWithID_DG(M1,outID)							#if each iteration of the small pole acheives a true flag, then the small list must be a subset of the large list
end

function getGainStats(M)
	av = 0
	min = M[1].gain
	max = M[1].gain
	for ID=1:length(M)
		av = av + M[ID].gain
		if M[ID].gain < min
			min = M[ID].gain
		end
		if M[ID].gain > max
			max = M[ID].gain
		end
	end
	av = av/length(M)
	return av, min, max
end

function getDeptStats(M)
	av = 0
	min = M[1].dept
	max = M[1].dept
	for ID=1:length(M)
		av = av + M[ID].dept
		if M[ID].dept < min
			min = M[ID].dept
		end
		if M[ID].dept > max
			max = M[ID].dept
		end
	end
	av = av/length(M)
	return av, min, max
end

function sizeStats(M)
	min = length(M[1].mem)
	max = length(M[1].mem)
	for ID=1:length(M)
		if length(M[ID].mem) < min
			min = length(M[ID].mem)
		end
		if length(M[ID].mem) > max
			max = length(M[ID].mem)
		end
	end
	return min,max
end

function findPolesWithID(M,IDs)
	outPoles = Any[]
	for i=1:length(IDs)
		outPoles = push!(outPoles,M[IDs[i]].mem)
	end
	return outPoles
end
function findPolesWithID_DG(M,IDs)
	outPoles = Any[]
	for i=1:length(IDs)
		outPoles = push!(outPoles,M[IDs[i]])
	end
	return outPoles
end
function supersetCheck(sP,bP)
	for bID=1:length(bP)
		bPm = sort(bP[bID].mem)
		for sID=1:length(sP)
			sPm = sort(sP[sID].mem)
			if issubset(sP,bP)
				return false
			end
		end
	end
	return true
end
function equals_PoleArray(P,A)
	return sort(A) == sort(P.mem)
end

function findArrayInPoles(A,M)
	A = sort(A)
	for i=1:length(M)
		if A == sort(M[i].mem)
			return true
		end
	end
	return false
end


function findArrayInClique(A,C)
	A = sort(A)
	for i=1:length(C)
		if A == sort(C[i])
			return true
		end
	end
	return false
end
function removeDupes(M)
	sort!(M, by=x->length(x.mem), rev=true)

    flag = falses(length(M))
	numpole = length(M)

    for x=numpole:-1:2
        for i=x-1:-1:1
            if issubset(M[x].mem, M[i].mem)
                flag[x] = true
                break
            end
        end
		# (numpole - x) % 10000 == 0 && (println("\t\tCompleted ", numpole - x, " / $numpole"))
    end

    deleteat!(M, findall(flag))

    return M
end

function isInOrIsSubsetIn(M,a)
	for i=1:length(M)
		pole = sort(M[i])
		a1 = sort(a)
		if issubset(a,pole)
			return true
		end
	end
	return false
end
function comparePolesEq(M1,M2)
	flag = true
	outID = []							#initalize the flag boolean
	for smallpoleID=1:length(M1)		#iterate through each member of the multipoles with less members
		flag = false					#set the flag to false before checking against the larger set
		sP = sort(M1[smallpoleID].mem)	#sort the current searched member of the smaller pole in number order. ie [3,2,1] becomes [1,2,3]
		for bigpoleID=1:length(M2)		#iterate through each member of the larger set for each member of the smaller set
			bP = sort(M2[bigpoleID].mem)#same sorting logic as the smaller poles for the current iteratation of the larger pole
			if sP == bP				#straight comparison of the two sorted arrays
				flag = true
				outID = append!(outID,smallpoleID)				#if the two arrays are equal, the smaller pole is a member of the larger poles
			end
		end
		if flag == false									#if a small pole array made it through the iteration of the larger poles without triggering the flag bool to true, the smaller collection
			outID = append!(outID,smallpoleID)				#	must not be a subset of the larger group, therefore return false.
		end												#	could be modified to return a list of non-subset poles and their dep/gain
	end
	return findPolesWithID(M1,outID)							#if each iteration of the small pole acheives a true flag, then the small list must be a subset of the large list
end
function printAdjOfPole(p, adj, c)
	num_var = size(c, 1)
	println(p.mem)
	println("Dep: ",p.dept)
	println("Gain:",p.gain)
	for i in p.mem
		for j in p.mem
			if c[i,j] <= 0.0
				print(adj[i,j])
				print(" ")
			else
				print(adj[i,j+num_var])
				print(" ")
			end
		end
		print("\n")
	end

	print("\n")
	for i in p.mem
		for j in p.mem
			print(c[i,j])
			print(" ")
		end
		print("\n")
	end

	print("\n\n\n")
end
function adjPrintAll(poles, adj, c)
	for p in poles
		printAdjOfPole(p,adj, c)
	end
end

function importAdj()
	outfile = "C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\New\\ArchiveOutput\\COMET_AdjM.csv"
	df = CSV.read(outfile, DataFrame,header = false)
    M = Matrix(df)
	c = importCSV()
	return M, c
end

function Lex(poles, sz)
	M1 = removeRedundancy(poles, sz)
	return M1
end



end
