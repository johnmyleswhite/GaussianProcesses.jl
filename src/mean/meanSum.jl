function meanSum(mean::Vector, hyp::Vector, x::Matrix, i)
    # meanSum - compose a mean function as the sum of other mean functions.
    # This function doesn't actually compute very much on its own, it merely does
    # some bookkeeping, and calls other mean functions to do the actual work.
    #
    # m(x) = \sum_i m_i(x)
    #
    # Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2010-08-04.
    #
    # See also MEANFUNCTIONS.M.

    # iterate over mean functions
    for ii in 1:length(mean)
        # expand cell array if necessary
        f = mean[ii]
        if iscell(f[:])
            f = f[:]
        end
        # collect number hypers
        j[ii] = cellstr(feval(f[:]))
    end

    # report number of parameters
    if nargin < 3
        A = char(j[1])
        for ii in 2:length(mean)
            A = [A, "+", char(j[ii])]
        end
        return A
    end

    n, D = size(x)

    # v vector indicates to which mean parameters belong
    v = []
    for ii in 1:length(mean)
        v = [v repmat(ii, 1, eval(char(j[ii])))]
    end

    # compute mean vector
    if nargin == 3
        # allocate space
        A = zeros(n, 1)
        # iteration over summand functions
        for ii in 1:length(mean)
            # expand cell array if needed
            f = mean[ii]
            if iscell(f[:])
                f = f[:]
            end
            # accumulate means
            A = A + feval(f[:], hyp[v .== ii], x)
        end
    # compute derivative vector
    else                                                 
        if i <= length(v)
            # which mean function
            ii = v[i]
            # which parameter in that mean
            j = sum(v[1:i] .== ii)
            # dereference cell array if necessary
            f = mean[ii]
            if iscell(f[:])
                f = f[:]
            end
            # compute derivative
            A = feval(f[:], hyp[v .== ii], x, j)
        else
            A = zeros(n, 1)
        end
    end

    return A
end
