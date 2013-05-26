function covScale(cov::Vector, hyp::Vector, x::Matrix, z::Matrix, i)
	# meanScale - compose a mean function as a scaled version of another one.
	#
	# k(x^p, x^q) = sf^2 * k_0(x^p, x^q)
	#
	# The hyperparameter is:
	#
	# hyp = [log(sf)]
	#
	# This function doesn't actually compute very much on its own, it merely does
	# some bookkeeping, and calls other mean function to do the actual work.
	#
	# Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2010-09-10.
	#
	# See also MEANFUNCTIONS.M.

	# report number of parameters
	if nargin < 3                                        
		K = [feval(cov[:]), "+1"]
		return K
	end

	# make sure, z exists
	if nargin < 4
		z = []
	end

	# signal variance
	n, D = size(x)
	sf2 = exp(2 * hyp[1])

	# covariances
	if nargin < 5
		K = sf2 * feval(cov[:], hyp[2:end], x, z)
	# derivatives
	else
		if i == 1
			K = 2 * sf2 * feval(cov[:], hyp[2:end], x, z)
	  	else
			K = sf2 * feval(cov[:], hyp[2:end], x, z, i - 1)
		end
	end

	return K
end
