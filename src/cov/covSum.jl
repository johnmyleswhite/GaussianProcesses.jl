function covSum(cov::Vector, hyp::Vector, x::Matrix, z::Matrix, i)
    # covSum - compose a covariance function as the sum of other covariance
    # functions. This function doesn't actually compute very much on its own, it
    # merely does some bookkeeping, and calls other covariance functions to do the
    # actual work.
    #
    # Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2010-09-10.
    #
    # See also COVFUNCTIONS.M.

    if length(cov) == 0
        error("We require at least one summand.")
    end

    # iterate over covariance functions
    for ii in 1:length(cov)
        # expand cell array if necessary
        f = cov[ii]
        if iscell(f[:])
            f = f[:]
        end
        # collect number hypers
        j[ii] = cellstr(feval(f[:]))
    end

    # report number of parameters
    if nargin < 3
        K = char(j[1])
        for ii in 2:length(cov)
            K = [K, "+", char(j[ii])]
        end
        return K
    end

    # make sure, z exists
    if nargin < 4
        z = []
    end

    n, D = size(x)

    # v vector indicates to which covariance parameters belong
    v = []

    for ii in 1:length(cov)
        v = [v repmat(ii, 1, eval(char(j(ii))))]
    end

    # covariances
    if nargin < 5
        K = 0.0
        # set default
        if nargin == 3
            z = []
        end
        # iteration over summand functions
        for ii in 1:length(cov)
            # expand cell array if necessary
            f = cov[ii]
            if iscell(f[:])
                f = f[:]
            end
            # accumulate covariances
            K = K + feval(f[:], hyp[v .== ii], x, z)
        end
    # derivatives
    else
        if i <= length(v)
            # which covariance function
            vi = v[i]
            # which parameter in that covariance
            j = sum(v[1:i] .== vi)
            f = cov[vi]
            # dereference cell array if necessary
            if iscell(f[:])
                f = f[:]
            end
            # compute derivative
            K = feval(f[:], hyp[v .== vi], x, z, j)
        else
            error("Unknown hyperparameter")
        end
    end

    return K
end
