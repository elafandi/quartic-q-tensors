function parse_distmesh(name)
	fp = open(name * "_p.txt", "r")
	V = countlines(fp)
	seek(fp, 0)
	P_raw = zeros(Float64, V, 2)
	for i=1:V
		P_raw[i,:] = map(k->parse(Float64,k), split(readline(fp)))
	end
	close(fp)
	
	tp = open(name * "_t.txt", "r")
	L = countlines(tp)
	seek(tp, 0)
	T_raw = zeros(Int, L, 3)
	for ell=1:L
		T_raw[ell,:] = map(k->Int(parse(Float64,k)), split(readline(tp)))
	end
	close(tp)
	
	edgeset = Set{Tuple{Int,Int}}()
	for ell=1:L
		a, b, c = T_raw[ell,:]
		for (p1, p2) in [(a,b), (a,c), (b,c)]
			if (p2, p1) in edgeset
				pop!(edgeset, (p2,p1))
			elseif (p1, p2) in edgeset
				pop!(edgeset, (p1,p2))
			else
				push!(edgeset, (p1,p2))
			end
		end
	end
	is_edgepoint = falses(V)
	for (p1, p2) in edgeset
		is_edgepoint[p1] = true
		is_edgepoint[p2] = true
	end
	
	N = V - count(is_edgepoint)
	remap_inv = [findall(s->!s, is_edgepoint); findall(is_edgepoint)]
	remap = zeros(Int, V)
	for i=1:V
		remap[remap_inv[i]] = i
	end
	P = zeros(Float64, V, 2)
	for i=1:V
		P[i,:] = P_raw[remap_inv[i],:]
	end
	T = map(i->remap[i], T_raw)
	
	return (P, T, N, V, L)
end