function meanScale(mean::Vector, hyp::Vector, x::Matrix, i)
	# meanScale - compose a mean function as a scaled version of another one.
	#
	# m(x) = a * m_0(x)
	#
	# The hyperparameter is:
	#
	# hyp = [ a ]
	#
	# This function doesn't actually compute very much on its own, it merely does
	# some bookkeeping, and calls other mean function to do the actual work.
	#
	# Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2010-07-15.
	#
	# See also MEANFUNCTIONS.M.

	# report number of parameters
	if nargin < 3
		A = [feval(mean[:]), "+1"]
		return A
	end

	n, D = size(x)
	a = hyp[1]
	# compute mean vector
	if nargin == 3
		A = a * feval(mean[:], hyp[2:end], x)
	else
		# compute derivative vector
		if i == 1
			A = feval(mean[:], hyp[2:end], x)
		else
			A = a * feval(mean[:], hyp[2:end], x, i - 1)
		end
	end

	return A
end
