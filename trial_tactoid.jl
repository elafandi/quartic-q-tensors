include("model2d.jl")
# include("triangulate_square.jl")
include("parse_distmesh.jl")

aBulk, bBulk, cBulk = -0.3, -4.0, 4.0

function init_fn(x,y)
	if x^2 + y^2 >= 0.3
		theta = atan(y, x)
		# n = [-sin(theta), cos(theta)]  # degree 1
		n = [-cos(theta), sin(theta)]  # degree -1
		# n = [1, 0]  # degree 0
		n ./= norm(n)
		sqrt(-2 * aBulk / cBulk) .* (n * n' - (1/2) .* I(2))
	else
		[0.0 0.0; 0.0 0.0]
	end
end

bound_fn(x,y) = init_fn(x,y)

# gridpoints = 20
# P, T, N, V, L = triangulate_square(2, 2, gridpoints)
# P, T, N, V, L = parse_distmesh("circle_362")
# P, T, N, V, L = parse_distmesh("circle_1452")
P, T, N, V, L = parse_distmesh("circle_5809")
timestep, max_time = 0.01, 90
params = ModelParams(
	P,
	T,
	N, V, L,
	max_time,
	timestep,
	1e-10,		# epsilon
	1e-1, 1e-3, 1e-3, 1e-3, 1e-3,		# L_1, L_2, L_3, L_4, L_5
	aBulk, bBulk, cBulk,
	1,		# M
	init_fn,
	bound_fn
)
fileprefix = "tactoid_results/degree_minus_one/tactoid_5809_-1"
Q = model(params, 2, false, fileprefix)