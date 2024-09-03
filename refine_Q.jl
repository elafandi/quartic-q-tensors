using LinearAlgebra
using Plots
include("triangulate_square.jl")

function refine_Q(fname_old, tris_per_side_old, fname_new, tris_per_side_new)
	Lx, Ly = 2, 2
	P_old, T_old, ins_old, N_old, L_old = triangulate_square(Lx, Ly, tris_per_side_old + 1)
	P_new, T_new, ins_new, N_new, L_new = triangulate_square(Lx, Ly, tris_per_side_new + 1)
	
	f = open(fname_old, "r")
	Q_old = map(ss->parse(Float64,strip(ss, ['[', ',', ']'])), split(readline(f)))
	close(f)
	Q_new = zeros(Float64, 2 * N_new)
	
	x_outside, y_outside = Float64[], Float64[]
	
	for i=1:N_new
		x,y = P_new[i,:]
		found = false
		for ell=1:L_old
			# is (x,y) inside Ω_ℓ?
			m1, n1 = P_old[T_old[ell,1],:]
			m2, n2 = P_old[T_old[ell,2],:]
			m3, n3 = P_old[T_old[ell,3],:]
			mnmat = [m1 m2 m3; n1 n2 n3; 1 1 1]
			bary = mnmat \ [x, y, 1]
			if minimum(bary) >= -1e-10
				# yes it is!
				found = true
				for alpha=1:2, k=1:3
					if T_old[ell,k] > N_old
						continue  # edge point; no Q
					end
					Q_new[2*(i-1)+alpha] += bary[k] * Q_old[2*(T_old[ell,k]-1)+alpha]
				end
				break
			end
		end
		if !found
			@show i, x, y
			push!(x_outside, x)
			push!(y_outside, y)
		end
	end
	
	open(fname_new, "w") do file
		write(file, string(Q_new))
	end
	
	scatter(x_outside, y_outside)
end