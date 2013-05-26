function meanMask(mean::Vector, hyp::Vector, x::Matrix, i)
	# Apply a mean function to a subset of the dimensions only. The subset can
	# either be specified by a 0/1 mask by a boolean mask or by an index set.
	#
	# This function doesn't actually compute very much on its own, it merely does
	# some bookkeeping, and calls another mean function to do the actual work.
	#
	# Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-01-16.
	#
	# See also MEANFUNCTIONS.M.

	# either a binary mask or an index set
	mask = fix(mean[1][:])

	# mean function to be masked
	mean = mean[2]

	# properly unwrap nested cell arrays
	if iscell(mean[:])
		mean = mean[:]
	end

	# number of hyperparameters of the full mean
	nh_string = feval(mean[:])

	# convert 1/0->index
	if max(mask) < 2 && length(mask) > 1
		mask = find(mask)
	end

	# masked dimension
	D = length(mask)

	# number of parameters
	if nargin < 3
		A = num2str(eval(nh_string))
		return
	end

	# check hyperparameters
	if eval(nh_string) != length(hyp)
		error("number of hyperparameters does not match size of masked data")
	end

	if nargin == 3
		A = feval(mean[:], hyp, x[:, mask])
	else
		# compute derivative
		A = feval(mean[:], hyp, x[:, mask], i)
	end

	return A
end
