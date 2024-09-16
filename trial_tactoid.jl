include("model2d.jl")
# include("triangulate_square.jl")
include("parse_distmesh.jl")

tact_degree = 0

function get_n(theta)
	if tact_degree == 1
		[-sin(theta), cos(theta)]
	elseif tact_degree == -1
		[-cos(theta), sin(theta)]
	else
		[1, 0]
	end
end

function init_fn(x,y)
	if x^2 + y^2 >= 0.3
		theta = atan(y, x)
		n = get_n(theta)
		n ./= norm(n)
		sqrt(-2 * aBulk / cBulk) .* (n * n' - (1/2) .* I(2))
	else
		[0.0 0.0; 0.0 0.0]
	end
end

bound_fn(x,y) = init_fn(x,y)

epsilon = 1e-10
L_1, L_2, L_3, L_4, L_5 = 1e-1, 1e-3, 1e-3, 1e-3, 1e-3
aBulk, bBulk, cBulk = -0.3, -4.0, 4.0

# gridpoints = 20
# P, T, N, V, L = triangulate_square(2, 2, gridpoints)
# P, T, N, V, L = parse_distmesh("circle_362")
# P, T, N, V, L = parse_distmesh("circle_1452")
P, T, N, V, L = parse_distmesh("circle_5809")

# If getting the energies, let max_time = 100 (it's enough)
timestep, max_time = 0.01, (tact_degree == 0 ? 40 : 1000)

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

# If getting the energies, comment the Q line and uncomment the ending lines.
if tact_degree == 1
	fileprefix = "tactoid_results/degree_one/tactoid_5809_1"
	Q = model(params, 3, false, fileprefix)
elseif tact_degree == -1
	fileprefix = "tactoid_results/degree_minus_one/tactoid_5809_-1"
	Q = model(params, 4, false, fileprefix)
elseif tact_degree == 0
	fileprefix = "tactoid_results/degree_zero/tactoid_5809_0"
	Q = model(params, 5, false, fileprefix)
end
# energies, Q = model(params, 2, true, fileprefix)
# open(fileprefix * "_energies.txt", "w") do file
	# write(file, string(energies))
# end
