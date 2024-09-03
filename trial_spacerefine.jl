include("model2d.jl")
include("triangulate_square.jl")

function init_fn(x, y)
	n = [
		x * (2-x) * y * (2-y),
		sin(pi*x) * sin(0.5*pi*y)
	]
	n * n' - (norm(n)^2 / 2) .* I(2)
end

function bound_fn(x, y)
	[0.0 0.0; 0.0 0.0]
end

max_time = 0.8
epsilon = 1e-10
L_1, L_2, L_3, L_4, L_5 = 1e-1, 1e-3, 1e-3, 1e-3, 1e-3
aBulk, bBulk, cBulk = -0.3, -4.0, 4.0

# for tris_per_side in [10, 20, 40, 80, 320]
for tris_per_side in [320]
	fileprefix = "space_refinement_results/" * string(tris_per_side) * "/"
	P, T, N, V, L = triangulate_square(2, 2, tris_per_side + 1)
	timestep = (tris_per_side == 320) ? max_time/25000 : max_time/1600
	params = ModelParams(
		P,
		T,
		N, V, L,
		max_time,
		timestep,
		epsilon,
		L_1, L_2, L_3, L_4, L_5,
		aBulk, bBulk, cBulk,
		1,		# M
		init_fn,
		bound_fn
	)
	if tris_per_side == 20
		energies, Q = model(params, 1, true, fileprefix)
		open(fileprefix * "energies.txt", "w") do file
			write(file, string(energies))
		end
	else
		Q = model(params, 1, false, fileprefix)
	end
end