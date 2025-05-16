using LinearAlgebra
using Plots
using SparseArrays
using ProgressMeter

struct ModelParams
	# Triangulation
	P::Matrix{Real}		# V×2; interior nodes are the top N×2
	T::Matrix{Int}		# L×3
	N::Int; V::Int; L::Int
	# Time domain
	max_time::Real
	timestep::Real
	epsilon::Real
	# Material constants
	L1::Real; L2::Real; L3::Real; L4::Real; L5::Real
	aBulk::Real; bBulk::Real; cBulk::Real
	M::Real
	# Initial and boundary conditions: takes in (x,y), returns 2x2 tensor
	init_fn::Function
	bound_fn::Function
end

function crossmat(v, mat)
	return [v[1] * mat[1,2] - v[2] * mat[1,1], v[1] * mat[2,2] - v[2] * mat[2,1]]
end

function model(params::ModelParams; get_energy=true, fileprefix=nothing, Q0=nothing, t0=0, graph_times::Vector{Float64}=Float64[], saveQ_times::Vector{Float64}=Float64[])
	# Step 1: Set up triangulation, global variables, and initial conditions

	# Triangulation
	P = params.P
	T = params.T
	N, V, L = params.N, params.V, params.L
	
	# HACKY MODIFICATION
	Pc, _, Nc, _, _ = parse_distmesh("circle_362")
	coarselist = Set{Int}()
	for i=1:Nc
		j = findmin(j->norm(Pc[i,:] - P[j,:]), 1:N)[2]
		push!(coarselist, j)
	end

	# Time domain
	max_time = params.max_time
	timestep = params.timestep
	epsilon = params.epsilon

	# Material constants
	L1, L2, L3, L4, L5 = params.L1, params.L2, params.L3, params.L4, params.L5
	aBulk, bBulk, cBulk = params.aBulk, params.bBulk, params.cBulk
	M = params.M
	s_0 = (bBulk + sqrt(bBulk^2 - 24 * aBulk * cBulk)) / (4 * cBulk)
	@assert 0 <= s_0 <= 1  # i guess it COULD be between -1/2 and 0 but that's not likely

	# Basis tensors of 2D symmetric traceless tensor space
	R = [[1 0; 0 -1], [0 1; 1 0]]
	
	# Boundary conditions
	Qbound = zeros(Float64, 2*(V-N))
	for i=1:V-N
		x,y = P[i+N,:]
		mat = params.bound_fn(x,y)
		Qbound[(i-1)*2 + 1] = mat[1,1]
		Qbound[(i-1)*2 + 2] = mat[1,2]
	end
	
	# Initial conditions
	if isnothing(Q0)
		Q0 = zeros(Float64, 2*N)
		for i=1:N
			x,y = P[i,:]
			mat = params.init_fn(x,y)
			Q0[(i-1)*2 + 1] = mat[1,1]
			Q0[(i-1)*2 + 2] = mat[1,2]
		end
		Q0 = [Q0; Qbound]
	end

	# Step 2: Find gradients of all psi's and integrals of all products of psi's

	mu = zeros(Float64, L, 3, 2)	# mu[ℓ,i,:] = gradient of ψ_(T[ℓ,i]) on Ω_ℓ
	Avec = zeros(Float64, L)  		# Avec[ℓ] = area of Ω_ℓ
	for ell = 1:L
		m1, n1 = P[T[ell,1],:]
		m2, n2 = P[T[ell,2],:]
		m3, n3 = P[T[ell,3],:]
		Avec[ell] = abs((m2 - m1) * (n3 - n1) - (m3 - m1) * (n2 - n1)) / 2
		mnmat = [m1 n1 1; m2 n2 1; m3 n3 1]
		for i=1:3
			mnvec = zeros(Float64, 3)
			mnvec[i] = 1
			a, b, _ = mnmat \ mnvec
			mu[ell,i,:] = [a, b]
		end
	end
	
	# Avec[ℓ] * K0 = area of Ω_ℓ
	# Avec[ℓ] * K1 = integral of ψ_(T[ℓ,i]) over Ω_ℓ
	# Avec[ℓ] * K2[i,j] = integral of ψ_(T[ℓ,i]) * ψ_(T[ℓ,j]) over Ω_ℓ
	# Avec[ℓ] * K4[i,j,k,z] = integral of product of four ψ's over Ω_ℓ
	K0 = 1
	K1 = 1 / 3
	K2 = (1/12) .* [2 1 1; 1 2 1; 1 1 2]
	# K3 = zeros(3,3,3)  # unnecessary because bBulk doesn't pop up in 2D
	# for i=1:3, j=1:3, k=1:3
		# if i == j == k
			# K3[i,j,k] = 1 / 10
		# elseif i == j || i == k || j == k
			# K3[i,j,k] = 1 / 30
		# else
			# K3[i,j,k] = 1 / 60
		# end
	# end
	K4 = zeros(3,3,3,3)
	for i=1:3, j=1:3, k=1:3, z=1:3
		if i == j == k == z
			K4[i,j,k,z] = 1 / 15
		elseif length(unique([i, j, k, z])) == 2
			if i == j == k || i == j == z || i == k == z || j == k == z
				K4[i,j,k,z] = 1 / 60
			else
				K4[i,j,k,z] = 1 / 90
			end
		else
			K4[i,j,k,z] = 1 / 180
		end
	end

	# Step 3: Find and store inner products of all pairs of internal basis functions

	A = spzeros(Float64, 2 * N, 2 * N)
	for ell = 1:L
		area = Avec[ell]
		for i_ind = 1:3, z_ind = 1:3
			i, z = T[ell,i_ind], T[ell,z_ind]
			if i > N || z > N
				continue
			end
			psiprod = area * K2[i_ind, z_ind]
			for alpha = 1:2
				y = alpha
				A[2*(z-1)+y,2*(i-1)+alpha] += 2 * psiprod
			end
		end
	end

	# Step 4: Set up function for fixed-point iteration
	
	# side1[ℓ,z_ind,β] represents R_β μ_{z,ℓ}
	# side2[ℓ,z_ind,α,β] represents R_α R_β μ_{z,ℓ}
	# sideC[ℓ,z_ind,β] represents μ_{z,ℓ} × R_β
	# each of these are 2-vecs, with components stored separately
	side11 = zeros(Float64, L, 3, 2)
	side12 = zeros(Float64, L, 3, 2)
	side21 = zeros(Float64, L, 3, 2, 2)
	side22 = zeros(Float64, L, 3, 2, 2)
	sideC1 = zeros(Float64, L, 3, 2)
	sideC2 = zeros(Float64, L, 3, 2)
	for ell = 1:L, z_ind=1:3, beta=1:2
		left = R[beta] * mu[ell,z_ind,:]
		side11[ell,z_ind,beta], side12[ell,z_ind,beta] = left
		for alpha=1:2
			right = R[alpha] * left
			side21[ell,z_ind,alpha,beta], side22[ell,z_ind,alpha,beta] = right
		end
		left = crossmat(mu[ell,z_ind,:], R[beta])
		sideC1[ell,z_ind,beta], sideC2[ell,z_ind,beta] = left
	end
	# mudot[ℓ,i_ind,j_ind] represents μ_{i,ℓ} • μ_{j,ℓ}
	mudot = zeros(Float64, L, 3, 3)
	for ell = 1:L, i_ind = 1:3, j_ind = 1:3
		mudot[ell,i_ind,j_ind] = dot(mu[ell,i_ind,:], mu[ell,j_ind,:])
	end

	function f(Q, X)
		H = zeros(Float64, 2 * N)
		
		@showprogress for ell = 1:L, y=1:2, z_ind=1:3
			area = Avec[ell]
			z = T[ell,z_ind]
			if z > N
				continue
			end
			
			# div and curl terms 1/6
			head_d1 = (L1 + 4*L3) * s_0^2/9
			head_c1 = (L2 + 4*L4) * s_0^2/9
			for j_ind=1:3
				j = T[ell,j_ind]
				for beta=1:2
					qjb = (Q[2*(j-1)+beta] + X[2*(j-1)+beta]) / 2
					# Div term
					dotted = side11[ell,j_ind,beta] * side11[ell,z_ind,y] +
							side12[ell,j_ind,beta] * side12[ell,z_ind,y]
					H[2*(z-1)+y] += head_d1 * qjb * dotted * area * K0
					# Curl term
					dotted = sideC1[ell,j_ind,beta] * sideC1[ell,z_ind,y] +
							sideC2[ell,j_ind,beta] * sideC2[ell,z_ind,y]
					H[2*(z-1)+y] += head_c1 * qjb * dotted * area * K0
				end
			end
			
			# All other curl terms will be zero, because we're multiplying an R tensor
			# basis function -- empty in the third row and column -- by the cross
			# product of two vectors in R2
			
			# term 2/6 and 3/6
			head_d2 = (L1 - 2*L3) * s_0/3
			head_d3 = (L1 - 2*L3) * s_0/3
			# head_c2 = (L2 - 2*L4) * s_0/3
			# head_c3 = (L2 - 2*L4) * s_0/3
			for j_ind=1:3, k_ind=1:3
				j, k = T[ell,j_ind], T[ell,k_ind]
				for beta=1:2, gamma=1:2
					qjb = (Q[2*(j-1)+beta] + X[2*(j-1)+beta]) / 2
					qkg = (Q[2*(k-1)+gamma] + X[2*(k-1)+gamma]) / 2
					# Div terms
					dotted = side11[ell,j_ind,beta] * side21[ell,z_ind,gamma,y] +
							side12[ell,j_ind,beta] * side22[ell,z_ind,gamma,y]
					H[2*(z-1)+y] += head_d2 * qjb * qkg * dotted * area * K1
					dotted = side11[ell,j_ind,beta] * side21[ell,k_ind,y,gamma] +
							side12[ell,j_ind,beta] * side22[ell,k_ind,y,gamma]
					H[2*(z-1)+y] += head_d3 * qjb * qkg * dotted * area * K1
				end
			end
			
			# term 4/6
			head_d4 = (L1 - 2*L3) * s_0/3
			# head_c4 = (L2 - 2*L4) * s_0/3
			for i_ind=1:3, j_ind=1:3
				i, j = T[ell,i_ind], T[ell,j_ind]
				for alpha=1:2, beta=1:2
					qiaqjb = (Q[2*(i-1)+alpha] * Q[2*(j-1)+beta] + 
								X[2*(i-1)+alpha] * X[2*(j-1)+beta]) / 2
					# Div term
					dotted = side21[ell,j_ind,alpha,beta] * side11[ell,z_ind,y] +
							side22[ell,j_ind,alpha,beta] * side12[ell,z_ind,y]
					H[2*(z-1)+y] += head_d4 * qiaqjb * dotted * area * K1
				end
			end
			
			# term 5/6 and 6/6
			head_d5 = L1 + L3
			head_d6 = L1 + L3
			# head_c5 = L2 + L4
			# head_c6 = L2 + L4
			for i_ind=1:3, j_ind=1:3, k_ind=1:3
				i, j, k = T[ell,i_ind], T[ell,j_ind], T[ell,k_ind]
				for alpha=1:2, beta=1:2, gamma=1:2
					qiaqjb = (Q[2*(i-1)+alpha] * Q[2*(j-1)+beta] + 
								X[2*(i-1)+alpha] * X[2*(j-1)+beta]) / 2
					qkg = (Q[2*(k-1)+gamma] + X[2*(k-1)+gamma]) / 2
					# Div terms
					dotted = side21[ell,j_ind,alpha,beta] * side21[ell,z_ind,gamma,y] +
							side22[ell,j_ind,alpha,beta] * side22[ell,z_ind,gamma,y]
					H[2*(z-1)+y] += head_d5 * qiaqjb * qkg * dotted * area * K2[i_ind,k_ind]
					dotted = side21[ell,j_ind,alpha,beta] * side21[ell,k_ind,y,gamma] +
							side22[ell,j_ind,alpha,beta] * side22[ell,k_ind,y,gamma]
					H[2*(z-1)+y] += head_d6 * qiaqjb * qkg * dotted * area * K2[i_ind,z_ind]
				end
			end
			
			# new energy term: |Q|^2 |∇Q|^2
			for i_ind=1:3, j_ind=1:3, k_ind=1:3
				i, j, k = T[ell,i_ind], T[ell,j_ind], T[ell,k_ind]
				for alpha=1:2
					beta = alpha
					gamma = y
					qiaqjb = (Q[2*(i-1)+alpha] * Q[2*(j-1)+beta] + 
								X[2*(i-1)+alpha] * X[2*(j-1)+beta]) / 2
					qkg = (Q[2*(k-1)+gamma] + X[2*(k-1)+gamma]) / 2
					H[2*(z-1)+y] += L5 * qiaqjb * qkg * 2 * 2 * (
							mudot[ell, k_ind, z_ind] * area * K2[i_ind, j_ind] +
							mudot[ell, i_ind, j_ind] * area * K2[k_ind, z_ind])
				end
			end
			
			# bulk energy, term a
			for i_ind=1:3
				i = T[ell,i_ind]
				alpha = y
				qia = (Q[2*(i-1)+alpha] + X[2*(i-1)+alpha]) / 2
				H[2*(z-1)+y] += 2 * aBulk * qia * 2 * area * K2[i_ind,z_ind]
			end
			
			# term c
			for i_ind=1:3, j_ind=1:3, k_ind=1:3
				i, j, k = T[ell,i_ind], T[ell,j_ind], T[ell,k_ind]
				for alpha=1:2
					beta = alpha
					gamma = y
					qiaqjb = (Q[2*(i-1)+alpha] * Q[2*(j-1)+beta] + 
								X[2*(i-1)+alpha] * X[2*(j-1)+beta]) / 2
					qkg = (Q[2*(k-1)+gamma] + X[2*(k-1)+gamma]) / 2
					H[2*(z-1)+y] += 2 * cBulk * qiaqjb * qkg * 2 *
									2 * area * K4[i_ind,j_ind,k_ind,z_ind]
				end
			end
		end
		[Q[1:2N] - M * timestep * (A \ H); Qbound]
	end
	
	function energy(Q)
		res = 0
		for ell = 1:L
			area = Avec[ell]
			# div and curl terms 1/3
			head_d = 1/2 * (L1 + 4*L3) * s_0^2/9
			head_c = 1/2 * (L2 + 4*L4) * s_0^2/9
			for j_ind=1:3, z_ind=1:3
				j, z = T[ell,j_ind], T[ell,z_ind]
				for beta=1:2, y=1:2
					qjb = Q[2*(j-1)+beta]
					qzy = Q[2*(z-1)+y]
					# Div term
					dotted = side11[ell,j_ind,beta] * side11[ell,z_ind,y] +
							side12[ell,j_ind,beta] * side12[ell,z_ind,y]
					res += head_d * qjb * qzy * dotted * area * K0
					# Curl term
					dotted = sideC1[ell,j_ind,beta] * sideC1[ell,z_ind,y] +
							sideC2[ell,j_ind,beta] * sideC2[ell,z_ind,y]
					res += head_c * qjb * qzy * dotted * area * K0
				end
			end
			
			# term 2/3
			head_d = (L1 - 2*L3) * s_0/3
			# head_c = (L2 - 2*L4) * s_0/3
			for j_ind=1:3, k_ind=1:3, z_ind=1:3
				j, k, z = T[ell,j_ind], T[ell,k_ind], T[ell,z_ind]
				for beta=1:2, gamma=1:2, y=1:2
					qjb = Q[2*(j-1)+beta]
					qkg = Q[2*(k-1)+gamma]
					qzy = Q[2*(z-1)+y]
					# Div terms
					dotted = side11[ell,j_ind,beta] * side21[ell,z_ind,gamma,y] +
							side12[ell,j_ind,beta] * side22[ell,z_ind,gamma,y]
					res += head_d * qjb * qkg * qzy * dotted * area * K1
				end
			end
			
			# term 3/3
			head_d = 1/2 * (L1 + L3)
			# head_c = 1/2 * (L2 + L4)
			for i_ind=1:3, j_ind=1:3, k_ind=1:3, z_ind=1:3
				i, j, k, z = T[ell,i_ind], T[ell,j_ind], T[ell,k_ind], T[ell,z_ind]
				for alpha=1:2, beta=1:2, gamma=1:2, y=1:2
					qia = Q[2*(i-1)+alpha]
					qjb = Q[2*(j-1)+beta]
					qkg = Q[2*(k-1)+gamma]
					qzy = Q[2*(z-1)+y]
					# Div terms
					dotted = side21[ell,j_ind,alpha,beta] * side21[ell,z_ind,gamma,y] +
							side22[ell,j_ind,alpha,beta] * side22[ell,z_ind,gamma,y]
					res += head_d * qia * qjb * qkg * qzy * dotted * area * K2[i_ind,k_ind]
				end
			end
			
			# new energy term: |Q|^2 |∇Q|^2
			for i_ind=1:3, j_ind=1:3, k_ind=1:3, z_ind=1:3
				i, j, k, z = T[ell,i_ind], T[ell,j_ind], T[ell,k_ind], T[ell,z_ind]
				for alpha=1:2, gamma=1:2
					beta = alpha
					y = gamma
					qia = Q[2*(i-1)+alpha]
					qjb = Q[2*(j-1)+beta]
					qkg = Q[2*(k-1)+gamma]
					qzy = Q[2*(z-1)+y]
					res += L5/2 * qia * qjb * qkg * qzy * 2 * 2 *
							mudot[ell, k_ind, z_ind] * area * K2[i_ind, j_ind]
				end
			end
			
			# bulk energy, term a
			for i_ind=1:3, z_ind=1:3
				i, z = T[ell,i_ind], T[ell,z_ind]
				for y=1:2
					alpha = y
					qia = Q[2*(i-1)+alpha]
					qzy = Q[2*(z-1)+y]
					res += aBulk * qia * qzy * 2 * area * K2[i_ind,z_ind]
				end
			end
			
			# term c
			for i_ind=1:3, j_ind=1:3, k_ind=1:3, z_ind=1:3
				i, j, k, z = T[ell,i_ind], T[ell,j_ind], T[ell,k_ind], T[ell,z_ind]
				for alpha=1:2, gamma=1:2
					beta = alpha
					y = gamma
					qia = Q[2*(i-1)+alpha]
					qjb = Q[2*(j-1)+beta]
					qkg = Q[2*(k-1)+gamma]
					qzy = Q[2*(z-1)+y]
					res += cBulk / 2 * qia * qjb * qkg * qzy * 2 *
									2 * area * K4[i_ind,j_ind,k_ind,z_ind]
				end
			end
		end
		res
	end
	
	function plot_eigenvec(Q, scale)
		plot(
			legend=false,
			xlim=[-1,1],
			ylim=[-1,1],
			aspect_ratio=:equal
		)
		for i=1:N
			mat = sum(alpha->R[alpha] .* Q[(i-1)*2+alpha], 1:2)
			if maximum(eigvals(mat)) == 0  # isotropic state
				continue
			end
			vind = findmax(eigvals(mat))[2]
			v = eigvecs(mat)[:,vind]
			v .*= scale / norm(v)
			plot!([P[i,1] - v[1] / 2, P[i,1] + v[1] / 2],
					[P[i,2] - v[2] / 2, P[i,2] + v[2] / 2], linecolor="blue")
		end
	end
	
	function plot_eigenval(Q)
		plot(
			xlim=[-1,1],
			ylim=[-1,1],
			aspect_ratio=:equal
		)
		for ell=1:L
			p1, p2, p3 = map(i->Tuple(P[i,:]), T[ell,:])
			Q1, Q2, Q3 = map(i->Q[2*i-1]*R[1] + Q[2*i]*R[2], T[ell,:])
			lambda = maximum(eigvals((Q1 + Q2 + Q3) ./ 3))
			tri = Shape([p1, p2, p3, p1])
			plot!(tri, fill_z = lambda, linealpha=0, label="")
		end
		for bound in [0.0, 0.2]  # the upper bound I saw in practice was around 0.2
			plot!(Shape(Tuple{Int64, Int64}[]), linealpha=0, label="", fill_z = bound)
		end
	end
	
	# HACKY MODIFICATION
	function plot_both(Q, scale, upper_bound)
		plot(
			# legend=false,
			xlim=[-1,1],
			ylim=[-1,1],
			aspect_ratio=:equal
		)
		@showprogress for ell=1:L
			p1, p2, p3 = map(i->Tuple(P[i,:]), T[ell,:])
			Q1, Q2, Q3 = map(i->Q[2*i-1]*R[1] + Q[2*i]*R[2], T[ell,:])
			lambda = maximum(eigvals((Q1 + Q2 + Q3) ./ 3))
			tri = Shape([p1, p2, p3, p1])
			plot!(tri, fill_z = lambda, linealpha=0, label="")
		end
		for bound in [0.0, upper_bound]
			plot!(Shape(Tuple{Int64, Int64}[]), linealpha=0, label="", fill_z = bound)
		end
		@showprogress for i in coarselist
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

	function Snorm(Q)
		res = 0.0
		for ell=1:L
			area = Avec[ell]
			for i_ind=1:3, j_ind=1:3
				i, j = T[ell,i_ind], T[ell,j_ind]
				for alpha=1:2
					beta = alpha
					res += 2 * Q[2*(i-1)+alpha] * Q[2*(j-1)+beta] * area * K2[i_ind,j_ind]
				end
			end
		end
		sqrt(res)
	end
	
	# Step 5: iterate!
	if get_energy
		energies = [energy(Q0)]
	end
	Q = Q0
	
	printing = false
	for current_time = t0+timestep:timestep:max_time
		oldX, newX = Q, f(Q, Q)
		iterations = 1
		curnorm = Snorm(oldX - newX)
		while curnorm > epsilon
			if iterations > 100
				error("no convergence in " * string(iterations) * " iterations")
			end
			if printing
				@show curnorm, iterations
			end
			oldX, newX = newX, f(Q, newX)
			iterations += 1
			curnorm = Snorm(oldX - newX)
		end
		if isnan(curnorm)
			error("curnorm is nan")
		end
		if printing
			@show curnorm, iterations
		end
		Q = newX
		
		if current_time in graph_times
			# HACKY MODIFICATION
			plot_both(Q, 0.075, 0.2)
			savefig(fileprefix * "_" * string(current_time) * ".png")
			# plot_eigenvec(Q, 0.025)
			# savefig(fileprefix * "_vec_" * string(current_time) * ".png")
			# plot_eigenval(Q)
			# savefig(fileprefix * "_val_" * string(current_time) * ".png")
		end
		if current_time in saveQ_times
			open(fileprefix * "_Q_" * string(current_time) * ".txt", "w") do file
				write(file, string(Q))  # to be turned into paper graphs later via make_tactoid_diagrams.jl
			end
		end
		
		if get_energy
			Qenergy = energy(Q)
			push!(energies, Qenergy)
		end
		
		if printing
			println("TIME:", current_time)
			println("ITERATIONS:", iterations)
			if get_energy
				println("||F(Q)|| = ", Qenergy)
			end
			println()
			# flush(stdout)
		else
			println("TIME:", current_time)
		end
	end
	if get_energy
		return (energies, Q)
	end
	return Q
end