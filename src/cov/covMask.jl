function covMask(cov::Vector, hyp::Vector, x::Matrix, z::Matrix, i)
    # Apply a covariance function to a subset of the dimensions only. The subset can
    # either be specified by a 0/1 mask by a boolean mask or by an index set.
    #
    # This function doesn't actually compute very much on its own, it merely does
    # some bookkeeping, and calls another covariance function to do the actual work.
    #
    # The function was suggested by Iain Murray, 2010-02-18 and is based on an
    # earlier implementation of his dating back to 2009-06-16.
    #
    # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2012-11-17.
    #
    # See also COVFUNCTIONS.M.

    # either a binary mask or an index set
    mask = fix(cov[1][:])

    # covariance function to be masked
    cov = cov[2]

    # properly unwrap nested cell arrays
    # TODO: Make this into a valid Julian type check
    if iscell(cov[:])
        cov = cov[:]
    end

    # number of hyperparameters of the full covariance
    nh_string = feval(cov[:])

    # convert 1/0 -> index
    if max(mask) < 2 && length(mask) > 1
        mask = find(mask)
    end

    # masked dimension
    D = length(mask)

    # number of parameters
    if nargin < 3
        K = num2str(eval(nh_string))
        return
    end

    # make sure, z exists
    if nargin < 4
        z = []
    end

    # determine mode
    xeqz = length(z) == 0
    dg = strcmp(z, "diag") && length(z) > 0

    if eval(nh_string) != length(hyp)
        # check hyperparameters
        error("number of hyperparameters does not match size of masked data")
    end

    # covariances
    if nargin < 5
        if dg
            K = feval(cov[:], hyp, x[:, mask], "diag")
        else
            if xeqz
                K = feval(cov[:], hyp, x[:, mask])
            else
                K = feval(cov[:], hyp, x[:, mask], z[:, mask])
            end
        end
    # derivatives
    else
        if i <= eval(nh_string)
            if dg
                K = feval(cov[:], hyp, x[:, mask], 'diag', i)
            else
                if xeqz
                    K = feval(cov[:], hyp, x[:, mask], [], i)
                else
                    K = feval(cov[:], hyp, x[:, mask], z[:, mask], i)
                end
            end
        else
            error('Unknown hyperparameter')
        end
    end

    return K
end
