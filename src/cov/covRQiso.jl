function covRQiso(hyp::Vector, x::Matrix, z::Matrix, i)
    # Rational Quadratic covariance function with isotropic distance measure. The
    # covariance function is parameterized as:
    #
    # k(x^p,x^q) = sf2 * [1 + (x^p - x^q)'*inv(P)*(x^p - x^q)/(2*alpha)]^(-alpha)
    #
    # where the P matrix is ell^2 times the unit matrix, sf2 is the signal
    # variance and alpha is the shape parameter for the RQ covariance. The
    # hyperparameters are:
    #
    # hyp = [ log(ell)
    #         log(sqrt(sf2))
    #         log(alpha) ]
    #
    # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
    #
    # See also COVFUNCTIONS.M.

    # report number of parameters
    if nargin < 2
        K = "3"
        return
    end

    # make sure, z exists
    if nargin < 3
        z = []
    end

    # determine mode
    xeqz = length(z) == 0
    dg = strcmp(z, "diag") && length(z) > 0

    ell = exp(hyp[1])
    sf2 = exp(2 * hyp[2])
    alpha = exp(hyp[3])

    # precompute squared distances
    if dg
        # vector kxx
        D2 = zeros(size(x, 1), 1)
    else
        if xeqz
            # symmetric matrix Kxx
            D2 = sq_dist(x' / ell)
        else
            # cross covariances Kxz
            D2 = sq_dist(x' / ell, z' / ell)
        end
    end

    # covariances
    if nargin < 4
        K = sf2 * ((1 + 0.5 * D2 / alpha).^(-alpha))
    # derivatives
    else
        # length scale parameter
        if i == 1
            K = sf2 * (1 + 0.5 * D2 / alpha).^(-alpha - 1) .* D2
        # magnitude parameter
        elseif i == 2
            K = 2 * sf2 * ((1 + 0.5 * D2 / alpha).^(-alpha))
        elseif i == 3
            K = (1 + 0.5 * D2 / alpha)
            K = sf2 * K.^(-alpha) .* (0.5 * D2 ./ K - alpha * log(K))
        else
            error("Unknown hyperparameter")
        end
    end

    return K
end
