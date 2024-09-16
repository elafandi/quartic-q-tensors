using LinearAlgebra
using SparseArrays
include("triangulate_square.jl")

function crossmat(v, mat)
	return [v[1] * mat[1,2] - v[2] * mat[1,1], v[1] * mat[2,2] - v[2] * mat[2,1]]
end

function get_energy(fname, gridpoints=400, fname2=nothing)
	Lx, Ly = 2, 2
	P, T, N, V, L = triangulate_square(Lx, Ly, gridpoints+1)
	
	s_0 = 1
	L1, L2, L3, L4, L5 = 1e-1, 1e-3, 1e-3, 1e-3, 1e-3
	aBulk, bBulk, cBulk = -0.3, -4.0, 4.0
	M = 1
	s_0 = (bBulk + sqrt(bBulk^2 - 24 * aBulk * cBulk)) / (4 * cBulk)
	
	R = [[1 0; 0 -1], [0 1; 1 0]]
	
	f = open(fname, "r")
	Q = map(ss->parse(Float64,strip(ss, ['[', ',', ']'])), split(readline(f)))
	close(f)
	
	if !isnothing(fname2)
		f2 = open(fname2, "r")
		Q2 = map(ss->parse(Float64,strip(ss, ['[', ',', ']'])), split(readline(f2)))
		close(f2)
		Q -= Q2
	end
	
	# Next part copied from model2d.jl
	
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
	
	if isnothing(fname2)
		energy(Q), Snorm(Q)
	else
		Snorm(Q)
	end
end