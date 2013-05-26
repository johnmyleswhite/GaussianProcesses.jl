function meanOne(hyp::Vector, x::Matrix, i)
	# One mean function. The mean function does not have any parameters.
	#
	# m(x) = 1
	#
	# Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-08-04.
	#
	# See also MEANFUNCTIONS.M.

	# report number of hyperparameters
	if nargin < 2
		A = "0"
		return
	end

	if nargin == 2
		# evaluate mean
		A = ones(size(x, 1), 1)
	else
		# derivative
		A = zeros(size(x, 1), 1)
	end

	return A
end
