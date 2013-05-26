function meanZero(hyp::Vector, x::Matrix, i)
	# Zero mean function. The mean function does not have any parameters.
	#
	# m(x) = 0
	#
	# Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-01-10.
	#
	# See also MEANFUNCTIONS.M.

	# report number of hyperparameters 
	if nargin < 2
		A = "0"
		return
	end

	# derivative and mean
	A = zeros(size(x, 1), 1)

	return A
end
