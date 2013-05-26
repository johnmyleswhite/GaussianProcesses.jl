function meanLinear(hyp::Vector, x::Matrix, i)
	# Linear mean function. The mean function is parameterized as:
	#
	# m(x) = sum_i a_i * x_i;
	#
	# The hyperparameter is:
	#
	# hyp = [ a_1
	#         a_2
	#         ..
	#         a_D ]
	#
	# Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-01-10.
	#
	# See also MEANFUNCTIONS.M.

	# report number of hyperparameters
	if nargin < 2
		A = "D"
		return
	end

	n, D = size(x)

	if any(size(hyp) .!= [D, 1])
		error("Exactly D hyperparameters needed.")
	end

	a = hyp

	if nargin == 2
		# evaluate mean
		A = x * a
	else
		if i <= D
			# derivative
	    	A = x[:, i]
		else
			A = zeros(n, 1)
		end
	end

	return A
end
