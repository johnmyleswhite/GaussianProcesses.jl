function covConst(hyp::Vector, x::Matrix, z::Matrix, i)
    # covariance function for a constant function.
    # The covariance function is parameterized as:
    #
    # k(x^p, x^q) = s2
    #
    # The scalar hyperparameter is:
    #
    # hyp = [log(sqrt(s2))]
    #
    # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
    #
    # See also COVFUNCTIONS.M.

    # report number of parameters
    if nargin < 2
        K = "1"
        return
    end

    # make sure, z exists
    if nargin < 3
        z = []
    end

    # determine mode
    xeqz = length(z) == 0
    dg = strcmp(z, "diag") && length(z) > 0

    # s2
    s2 = exp(2 * hyp)
    n = size(x, 1)

    # vector kxx
    if dg
        K = s2 * ones(n, 1)
    else
        # symmetric matrix Kxx
        if xeqz
            K = s2 * ones(n)
        # cross covariances Kxz
        else
            K = s2 * ones(n, size(z, 1))
        end
    end

    # derivatives
    if nargin > 3
        if i == 1
            K = 2 * K
        else
            error("Unknown hyperparameter")
        end
    end

    return K
end
