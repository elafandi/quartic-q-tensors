include("parse_distmesh.jl")
using LinearAlgebra
using Images
using Plots

fileprefix = "tactoid_results/degree_one/tactoid_5809_1"
times = [0.01, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 180.0, 270.0]

# fileprefix = "tactoid_results/degree_minus_one/tactoid_5809_-1"
# times = [0.01, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 60.0, 90.0]

# fileprefix = "tactoid_results/degree_zero/tactoid_5809_0"
# times = [0.01, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]

P, T, N, V, L = parse_distmesh("circle_5809")
R = [[1 0; 0 -1], [0 1; 1 0]]
Pc, _, Nc, _, _ = parse_distmesh("circle_362")
coarselist = Set{Int}()
for i=1:Nc
	j = findmin(j->norm(Pc[i,:] - P[j,:]), 1:N)[2]
	push!(coarselist, j)
end

function plot_both(Q, scale, upper_bound)
	plot(
		xlim=[-1,1],
		ylim=[-1,1],
		# This one is going in the paper -- it has to be smaller
		size=(265,250), legend=false
		# aspect_ratio=:equal
	)
	for ell=1:L
		p1, p2, p3 = map(i->Tuple(P[i,:]), T[ell,:])
		Q1, Q2, Q3 = map(i->Q[2*i-1]*R[1] + Q[2*i]*R[2], T[ell,:])
		lambda = maximum(eigvals((Q1 + Q2 + Q3) ./ 3))
		tri = Shape([p1, p2, p3, p1])
		plot!(tri, fill_z = lambda, linealpha=0, label="")
	end
	for bound in [0.0, upper_bound]
		plot!(Shape(Tuple{Int64, Int64}[]), linealpha=0, label="", fill_z = bound)
	end
	for i in coarselist  # and because it's smaller, we can only plot the coarse vectors
		mat = sum(alpha->R[alpha] .* Q[(i-1)*2+alpha], 1:2)
		if maximum(eigvals(mat)) == 0  # isotropic state
			continue
		end
		vind = findmax(eigvals(mat))[2]
		v = eigvecs(mat)[:,vind]
		v .*= scale / norm(v)
		plot!([P[i,1] - v[1] / 2, P[i,1] + v[1] / 2],
				[P[i,2] - v[2] / 2, P[i,2] + v[2] / 2], linecolor="blue", label="")
	end
end

Qvec = Vector{Float64}[]
# upper_bound = 0.0
for t in times
	open(fileprefix * "_Q_" * string(t) * ".txt") do file
		w = readline(file)
		Q = map(ss->parse(Float64,strip(ss, ['[', ',', ']'])), split(w))
		push!(Qvec, Q)
		# for ell=1:L
			# Q1, Q2, Q3 = map(i->Q[2*i-1]*R[1] + Q[2*i]*R[2], T[ell,:])
			# lambda = maximum(eigvals((Q1 + Q2 + Q3) ./ 3))
			# global upper_bound = max(upper_bound, lambda)
		# end
	end
end
# @show upper_bound
upper_bound = 0.2  # better to have it consistent; it's just below for each anyway
for i=1:9
	Q, t = Qvec[i], times[i]
	plot_both(Q, 0.05, upper_bound)
	fname = fileprefix * "_both_" * string(t) * ".png"
	savefig(fname)
	A = load(fname)  # there's probably a more efficient way but meh
	A = A[1:end-15, 10:end]
	save(fname, A)
end
