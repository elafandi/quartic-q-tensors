include("model2d.jl")
include("triangulate_square.jl")

function init_fn(x, y)
	n = [
		x * (2-x) * y * (2-y),
		sin(pi*x) * sin(0.5*pi*y)
	]
	n * n' - (norm(n)^2 / 2) .* I(2)
end

bound_fn(x, y) = Float64[0 0; 0 0]

epsilon = 1e-10
L_1, L_2, L_3, L_4, L_5 = 1e-1, 1e-3, 1e-3, 1e-3, 1e-3
aBulk, bBulk, cBulk = -0.3, -4.0, 4.0

function is_convergent(tris_per_side, dt, max_time)
	P, T, N, V, L = triangulate_square(2, 2, tris_per_side + 1)
	params = ModelParams(
		P,
		T,
		N, V, L,
		max_time,
		dt,
		epsilon,
		L_1, L_2, L_3, L_4, L_5,
		aBulk, bBulk, cBulk,
		1,		# M
		init_fn,
		bound_fn
	)
	try
		model(params, 2, false)
	catch
		return false
	end
	return true
end

dt_L = 0.1
threshold = 1e-6

for tris_per_side in [10, 20, 40, 80, 160, 320]
	global dt_L
	timesteps_needed = (tris_per_side == 320 ? 100 : 50)
	while !is_convergent(tris_per_side, dt_L, dt_L * timesteps_needed)
		dt_L /= 2
	end
	dt_R = dt_L * 2
	while log(dt_R) - log(dt_L) > threshold
		dt_M = (dt_L + dt_R) / 2
		@show tris_per_side, dt_L, dt_R, log(dt_R) - log(dt_L)
		if is_convergent(tris_per_side, dt_M, dt_M * timesteps_needed)
			dt_L = dt_M
		else
			dt_R = dt_M
		end
	end
	@show tris_per_side, dt_L, dt_R
	@assert is_convergent(tris_per_side, dt_L, dt_L * timesteps_needed * 2)
	open("cfl_results.txt", "a") do file
		s = string(tris_per_side) * ": " * string(dt_L) * " " * string(dt_R) * "\n"
		print(s)
		write(file, s)
	end
end