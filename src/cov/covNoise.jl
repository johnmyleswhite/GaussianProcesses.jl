function covNoise(hyp::Vector, x::Matrix, z::Matrix, i)
    # Independent covariance function, ie "white noise",
    # with specified variance. The covariance function is specified as:
    #
    # k(x^p,x^q) = s2 * \delta(p, q)
    #
    # where s2 is the noise variance and \delta(p,q) is a Kronecker delta
    # function which is 1 iff p=q and zero otherwise. Two data points
    # p and q are considered equal if their norm is less than 1e-9.
    # The hyperparameter is
    #
    # hyp = [log(sqrt(s2))]
    #
    # For more help on design of covariance functions, try "help covFunctions".
    #
    # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2012-12-19.
    #
    # See also COVFUNCTIONS.M.#

    # threshold on the norm when two vectors are considered to be equal
    tol = 1e-9

    # report number of parameters
    if nargin < 2
        K = "1"
        return K
    end

    # make sure, z exists
    if nargin < 3
        z = []
    end

    # determine mode
    xeqz = length(z) == 0
    dg = strcmp(z, "diag") && length(z) > 0

    if ndims(x) == ndims(z) && all(size(x) == size(z))
        xeqz = norm(x - z, Inf) < tol
    end

    # noise variance
    n = size(x, 1)
    s2 = exp(2 * hyp)

    # precompute raw
    if dg
        # vector kxx
        K = ones(n, 1)
    else
        if xeqz
            # symmetric matrix Kxx
            K = eye(n)
        else
            # cross covariances Kxz
            K = double(sq_dist(x', z') < tol * tol)
        end
    end

    if nargin < 4
        # covariances
        K = s2 * K
    else
        # derivatives
        if i == 1
            K = 2 * s2 * K
        else
            error("Unknown hyperparameter")
        end
    end

    return K
end
