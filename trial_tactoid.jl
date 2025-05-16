include("model2d.jl")
# include("triangulate_square.jl")
include("parse_distmesh.jl")

tact_degree = 2
max_time = 270

timestep = 0.01
epsilon = 1e-10

L_1, L_2, L_3, L_4, L_5 = 1e-1, 1e-3, 1e-3, 1e-3, 1e-3
aBulk, bBulk, cBulk = -0.3, -4.0, 4.0
M = 1

if tact_degree == 1
	dir = "tactoid_results/degree_one_L1=" * string(L_1) * "/"
	fileprefix = dir * "tactoid_5809_1"
	saveQ_times = [timestep; collect(5:5:max_time)]
	graph_times = [timestep; collect(5:5:max_time)]
elseif tact_degree == -1
	dir = "tactoid_results/degree_minus_one_L1=" * string(L_1) * "/"
	fileprefix = dir * "tactoid_5809_-1"
	saveQ_times = [timestep; collect(5:5:max_time)]
	graph_times = [timestep; collect(5:5:max_time)]
elseif tact_degree == 0
	dir = "tactoid_results/degree_zero_L1=" * string(L_1) * "/"
	fileprefix = dir * "tactoid_5809_0"
	saveQ_times = [timestep; collect(5:5:max_time)]
	graph_times = [timestep; collect(5:5:max_time)]
elseif tact_degree == 2  # technically incorrect notation but this makes it easier
	dir = "tactoid_results/reversed_L1=" * string(L_1) * "/"
	fileprefix = dir * "tactoid_5809_R"
	saveQ_times = [timestep; collect(5:5:max_time)]
	graph_times = [timestep; collect(5:5:max_time)]
end

# file = open(fileprefix * "_Q_270.0.txt")
# w = readline(file)
# Q0 = map(ss->parse(Float64,strip(ss, ['[', ',', ']'])), split(w))

mkpath(dir)

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
	if (x^2 + y^2 >= 0.3 && tact_degree in [-1, 0, 1]) || (x^2 + y^2 <= 0.3 && tact_degree == 2)
		theta = atan(y, x)
		n = get_n(theta)
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

params = ModelParams(
	P,
	T,
	N, V, L,
	max_time,
	timestep,
	epsilon,
	L_1, L_2, L_3, L_4, L_5,
	aBulk, bBulk, cBulk,
	M,
	init_fn,
	bound_fn
)

# leftres, rightres = model(params)
energies, Q = model(params, get_energy=true, fileprefix=fileprefix, saveQ_times=saveQ_times, graph_times=graph_times)
open(dir * "energies_0_to_270.txt", "w") do file
	print(file, energies)
end