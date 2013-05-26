function covSEard(hyp::Vector, x::Matrix, z::Matrix, i)
    # Squared Exponential covariance function with Automatic Relevance Detemination
    # (ARD) distance measure. The covariance function is parameterized as:
    #
    # k(x^p, x^q) = sf2 * exp(-(x^p - x^q)' * inv(P) * (x^p - x^q) / 2)
    #
    # where the P matrix is diagonal with ARD parameters ell_1^2,...,ell_D^2, where
    # D is the dimension of the input space and sf2 is the signal variance. The
    # hyperparameters are:
    #
    # hyp = [ log(ell_1)
    #         log(ell_2)
    #          .
    #         log(ell_D)
    #         log(sqrt(sf2)) ]
    #
    # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
    #
    # See also COVFUNCTIONS.M.

    # report number of parameters
    if nargin < 2
        K = "(D+1)"
        return
    end

    # make sure, z exists
    if nargin < 3
        z = []
    end

    # determine mode
    xeqz = length(z) == 0
    dg = strcmp(z, "diag") && length(z) > 0

    n, D = size(x)
    # characteristic length scale
    ell = exp(hyp[1:D])
    # signal variance
    sf2 = exp(2 * hyp[D + 1])

    # precompute squared distances
    if dg
        # vector kxx
        K = zeros(size(x, 1), 1)
    else
        if xeqz
            # symmetric matrix Kxx
            K = sq_dist(diag(1 ./ ell) * x')
        else
            # cross covariances Kxz
            K = sq_dist(diag(1 ./ ell) * x', diag(1 ./ ell) * z')
        end
    end

    # covariance
    K = sf2 * exp(-K / 2)

    # derivatives
    if nargin > 3
        # length scale parameters
        if i <= D
            if dg
                K = K * 0.0
            else
                if xeqz
                    K = K .* sq_dist(x[:, i]' / ell[i])
                else
                    K = K .* sq_dist(x[:, i]' / ell[i], z[:, i]' / ell[i])
                end
            end
        # magnitude parameter
        elseif i == D + 1
            K = 2 * K
        else
            error("Unknown hyperparameter")
        end
    end

    return K
end
