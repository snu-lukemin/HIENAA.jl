"""
leastsquares interpolates d+2 points to a polynomial of degree d so that they have identical size of errors using least squares. 
"""
@views function leastsquares(f::Function, xvals::Vector{BigFloat}, a::Real, b::Real)::Vector{BigFloat}
	d = length(xvals)
	matrix = Array{BigFloat}(undef, d, d)

	vector = @. 2f((b - a) / 2 * (xvals + 1) + a) / (b - a)
	@. matrix[:, 1] = 1
	@. matrix[:, 2] = xvals

	for i ∈ 3:d-1
		@. matrix[:, i] = 2 * xvals * matrix[:, i-1] - matrix[:, i-2]
	end

	@. matrix[:, d] = (-1)^(1:d)

	matrix \ vector
end

"""
findnode finds the piecewise chebyshev node which minimises the error.
This is a Julia implementation of the node finding algorithm by Han et al.(https://eprint.iacr.org/2019/688)
"""
function findnode(interval::Vector{NTuple{2, BigFloat}}, n::Int, a::BigFloat, b::BigFloat)::Vector{BigFloat}
	m = length(interval)
	grid = Vector{BigFloat}(undef, n)

	if n < m
		jump = BigFloat(0.0)
		for i ∈ eachindex(interval)
			jump += interval[i][2] - interval[i][1]
		end

		deg = zeros(Int64, m)

		j = 1
		diff = interval[1][2] - interval[1][1]
		val = 0
		for i ∈ 1:n
			val += jump
			while diff < jump
				val -= interval[j][2] - interval[j][1]
				j += 1
				diff += interval[j][2] - interval[j][1]
			end
			deg[j] = 1
			grid[i] = interval[j][1] + val
			diff -= jump
		end
	else
		deg = optinode(interval, n)
		cnt = 0
		for i ∈ 1:m
			len = deg[i]
			centre = (interval[i][2] + interval[i][1]) / 2
			halfwidth = (interval[i][2] - interval[i][1]) / 2
			for i ∈ 1:len
				grid[cnt+i] = centre + halfwidth * cospi((2(len - i + 1) - 1) / (2 * len))
			end

			cnt += len
		end
	end

    for i = eachindex(grid)
        grid[i] = 2(grid[i] - a) / (b - a) - 1
    end

    grid
end

"""
optinode returns the optimal number of nodes for each interval.
Currently, it can only compute the case where the length of intervals are the same.
"""
function optinode(interval::Vector{NTuple{2, BigFloat}}, n::Int)::Vector{Int64}
	m = length(interval)
	centres = Vector{BigFloat}(undef, m)
	mininterval = Inf
	for i ∈ 1:m
		centres[i] = (interval[i][2] + interval[i][1]) / 2
		if mininterval > interval[i][2] - interval[i][1]
			mininterval = interval[i][2] - interval[i][1]
		end
	end
	dev = -floor(Int64, log2(mininterval))
	tot_deg = m
    deg = ones(Int64, m)

	err = 1 / big(1 << dev)
    tmp = m * log2(π) + log2(err)
    for i = 1:m
        tmp -= log2(i)
    end
    bdd = fill(tmp, m)
    for i = 1:m
        for j = 1:i-1
            bdd[i] += log2(-centres[j] + (centres[i] + err))
        end
        for j = i+1:m
            bdd[i] += log2(centres[j] - (centres[i] + err))
        end
    end

	while tot_deg < n
		maxi = findfirst(x -> x == maximum(bdd), bdd)
        for i = 1:m
            bdd[i] += log2(abs(centres[i] - centres[maxi]) + err) + log2(2π) - log2(tot_deg + 1)
        end
		bdd[maxi] -= 1
		tot_deg += 1
		deg[maxi] += 1
	end

	deg
end

function approximate(f::Function, interval::Vector{NTuple{2, BigFloat}}, d::Int64, a::Real, b::Real)::Vector{BigFloat}
	setprecision(256)

    a, b = BigFloat(a), BigFloat(b)
	node = findnode(interval, d + 2, a, b)    
	res = leastsquares(f, node, a, b)
	pop!(res)

	res
end