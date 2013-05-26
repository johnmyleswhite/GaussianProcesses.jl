function meanPow(mean::Vector, hyp::Vector, x::Matrix, i)
	# meanPow - compose a mean function as the power of another mean function.
	#
	# The degree d has to be a strictly positive integer.
	#
	# m(x) = m_0(x) ^ d
	#
	# This function doesn't actually compute very much on its own, it merely does
	# some bookkeeping, and calls other mean function to do the actual work.
	#
	# Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2010-06-18.
	#
	# See also MEANFUNCTIONS.M.

	# positive integer degree
	d = mean[1]
	d = abs(floor(d))
	d = max(d, 1)

	mean = mean[2]
	if !iscell(mean)
		mean = {mean}
	end

	# report number of parameters
	if nargin < 3
		A = feval(mean[:])
		return A
	end

	n, D = size(x)
	# compute mean vector
	if nargin == 3
		A = feval(mean[:], hyp, x).^d
	# compute derivative vector
	else
		A = (d * feval(mean[:], hyp, x).^(d - 1)) .* feval(mean[:], hyp, x, i)
	end

	return A
end
