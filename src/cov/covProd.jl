function covProd(cov, hyp::Vector, x::Matrix, z::Matrix, i)
    # covProd - compose a covariance function as the product of
    # other covariance functions. This function doesn't actually
    # compute very much on its own, it merely does some bookkeeping,
    # and calls other covariance functions to do the actual work.
    #
    # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
    #
    # See also COVFUNCTIONS.M.

    if length(cov) == 0
        error("We require at least one factor.")
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
        v = [v repmat(ii, 1, eval(char(j[ii])))]
    end

    # covariances
    if nargin < 5
        K = 1
        # set default
        if nargin == 3
            z = []
        end
        # iteration over factor functions
        for ii in 1:length(cov)
            # expand cell array if necessary
            f = cov[ii]
            if iscell(f[:])
                f = f[:]
            end
            # accumulate covariances
            K = K .* feval(f[:], hyp(v .== ii), x, z)
        end
    # derivatives
    else
        if i <= length(v)
            K = 1
            # which covariance function
            vi = v[i]
            # which parameter in that covariance
            j = sum(v[1:i] .== vi)
            #  iteration over factor functions
            for ii in 1:length(cov)
                % expand cell if necessary
                f = cov[ii]
                if iscell(f[:])
                    f = f[:]
                end
                if ii == vi
                    # accumulate covariances
                    K = K .* feval(f[:], hyp[v .== ii], x, z, j)
                else
                    # accumulate covariances
                    K = K .* feval(f[:], hyp[v .== ii], x, z)
                end
            end
        else
            error("Unknown hyperparameter")
        end
    end

    return K
end
