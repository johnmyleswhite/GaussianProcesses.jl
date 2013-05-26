function covLINone(hyp::Vector, x::Matrix, z::Matrix, i)
    # Linear covariance function with a single hyperparameter.
    # The covariance function is parameterized as:
    #
    # k(x^p, x^q) = (x^p' * x^q + 1) / t2;
    #
    # where the P matrix is t2 times the unit matrix. The second term plays the
    # role of the bias. The hyperparameter is:
    #
    # hyp = [log(sqrt(t2))]
    #
    # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
    #
    # See also COVFUNCTIONS.M.#

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

    # t2 inverse
    it2 = exp(-2 * hyp)

    # precompute inner products
    if dg
        # vector kxx
        K = sum(x .* x, 2)
    else
        if xeqz
            # symmetric matrix Kxx
            K = x * x'
        else
            # cross covariances Kxz
            K = x * z'
        end
    end

    if nargin < 4
        # covariances
        K = it2 * (1 + K)
    else
        # derivatives
        if i == 1
            K = -2 * it2 * (1 + K)
        else
            error("Unknown hyperparameter")
        end
    end

    return K
end
