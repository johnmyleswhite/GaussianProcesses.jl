function meanConst(hyp::Vector, x::Matrix, i)
	# Constant mean function. The mean function is parameterized as:
	#
	# m(x) = c
	#
	# The hyperparameter is:
	#
	# hyp = [c]
	#
	# Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-08-04.
	#
	# See also MEANFUNCTIONS.M.

	# report number of hyperparameters 
	if nargin < 2
		A = "1"
		return
	end

	if length(hyp) != 1
		error("Exactly one hyperparameter needed.")
	end

	c = hyp
	if nargin == 2
		# evaluate mean
		A = c * ones(size(x, 1), 1)
	else
		if i == 1
			# derivative
			A = ones(size(x, 1), 1)
		else
			A = zeros(size(x, 1), 1)
		end
	end

	return A
end
