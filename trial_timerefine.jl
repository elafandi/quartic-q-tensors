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

tris_per_side = 30
max_time = 0.8
epsilon = 1e-10
L_1, L_2, L_3, L_4, L_5 = 1e-1, 1e-3, 1e-3, 1e-3, 1e-3
aBulk, bBulk, cBulk = -0.3, -4.0, 4.0

for t_count in [200, 400, 800, 1600, 3200, 80000]
	timestep = 0.8 / t_count
	fileprefix = "time_refinement_results/" * string(t_count) * "/"
	P, T, N, V, L = triangulate_square(2, 2, tris_per_side + 1)
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
	Q = model(params, 1, false, fileprefix)
end