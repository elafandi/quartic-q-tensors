function triangulate_square(Lx, Ly, n)
	w = n - 2
	N = w^2  # number of interior nodes
	V = w^2 + 4*w + 2  # number of nodes (including boundary)
	L = 2 * (w+1)^2 - 2  # number of triangles

	P = zeros(Float64, V, 2)  # coordinates of points
	for y=1:w, x=1:w
		P[(y-1)*w + x, :] = [x * Lx / (w+1), y * Ly / (w+1)]
	end
	for t=1:w
		P[w^2 + t, :] = [t * Lx / (w+1), 0]
		P[w^2 + w + t, :] = [Lx, t * Ly / (w+1)]
		P[w^2 + 2 * w + t, :] = [t * Lx / (w+1), Ly]
		P[w^2 + 3 * w + t, :] = [0, t * Ly / (w+1)]
	end
	P[end-1,:] = [Lx,Ly]
	P[end,:] = [0,0]

	T = zeros(Int, L, 3)  # indices of points forming each triangle
	tind = 1
	for y=1:w-1, x=1:w-1
		i = (y-1) * w + x
		T[tind,:] = [i, i + w + 1, i + w]
		T[tind+1,:] = [i, i + 1, i + w + 1]
		tind += 2
	end
	for x=1:w-1
		T[tind,:] = [x+1, x, w^2+x]
		T[tind+1,:] = [x+1, w^2+x, w^2+x+1]
		tind += 2
	end
	for y=1:w-1
		T[tind,:] = [y*w, w^2+w+y, w^2+w+y+1]
		T[tind+1,:] = [(y+1)*w, y*w, w^2+w+y+1]
		tind += 2
	end
	for x=1:w-1
		T[tind,:] = [w^2-w+x, w^2+2*w+x+1, w^2+2*w+x]
		T[tind+1,:] = [w^2-w+x, w^2-w+x+1, w^2+2*w+x+1]
		tind += 2
	end
	for y=1:w-1
		T[tind,:] = [(y-1)*w+1, y*w+1, w^2+3*w+y]
		T[tind+1,:] = [y*w+1, w^2+3*w+y+1, w^2+3*w+y]
		tind += 2
	end
	T[end-5,:] = [w, w^2+w, w^2+w+1]
	T[end-4,:] = [w^2, w^2+2*w, w^2+4*w+1]
	T[end-3,:] = [w^2, w^2+4*w+1, w^2+3*w]
	T[end-2,:] = [w^2-w+1, w^2+2*w+1, w^2+4*w]
	T[end-1,:] = [1, w^2+3*w+1, w^2+4*w+2]
	T[end,:] = [1, w^2+4*w+2, w^2+1]
	return (P, T, N, V, L)
end