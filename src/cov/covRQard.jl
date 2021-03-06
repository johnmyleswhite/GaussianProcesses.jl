function covRQard(hyp::Vector, x::Matrix, z::Matrix, i)
    # Rational Quadratic covariance function with Automatic Relevance Determination
    # (ARD) distance measure. The covariance function is parameterized as:
    #
    # k(x^p,x^q) = sf2 * [1 + (x^p - x^q)'*inv(P)*(x^p - x^q)/(2*alpha)]^(-alpha)
    #
    # where the P matrix is diagonal with ARD parameters ell_1^2,...,ell_D^2, where
    # D is the dimension of the input space, sf2 is the signal variance and alpha
    # is the shape parameter for the RQ covariance. The hyperparameters are:
    #
    # hyp = [ log(ell_1)
    #         log(ell_2)
    #          ..
    #         log(ell_D)
    #         log(sqrt(sf2))
    #         log(alpha) ]
    #
    # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-08-04.
    #
    # See also COVFUNCTIONS.M.

    # report number of parameters
    if nargin < 2
        K = "(D + 2)"
        return K
    end

    # make sure, z exists
    if nargin < 3
        z = []
    end

    # determine mode
    xeqz = length(z) == 0
    dg = strcmp(z, "diag") && length(z) > 0

    n, D = size(x)
    ell = exp(hyp[1:D])
    sf2 = exp(2 * hyp[D + 1])
    alpha = exp(hyp[D + 2])

    # precompute squared distances
    if dg
        # vector kxx
        D2 = zeros(size(x, 1), 1)
    else
        if xeqz
            # symmetric matrix Kxx
            D2 = sq_dist(diag(1 ./ ell) * x')
        else
            # cross covariances Kxz
            D2 = sq_dist(diag(1 ./ ell) * x', diag(1 ./ ell) * z')
        end
    end

    # covariances
    if nargin < 4
        K = sf2 * ((1 + 0.5 * D2 / alpha).^(-alpha))
    # derivatives
    else                                                               
        # length scale parameter
        if i <= D
            if dg
                K = D2 * 0.0
            else
                if xeqz
                    K = sf2 * (1 + 0.5 * D2 / alpha).^(-alpha - 1) .* sq_dist(x[:, i]' / ell[i])
                else
                    K = sf2 * (1 + 0.5 * D2 / alpha).^(-alpha - 1) .* sq_dist(x[:, i]' / ell[i], z[:,i]' / ell[i])
                end
            end
        # magnitude parameter
        elseif i == D + 1
            K = 2 * sf2 * ((1 + 0.5 * D2 / alpha).^(-alpha))
        elseif i == D + 2
            K = (1 + 0.5 * D2 / alpha)
            K = sf2 * K.^(-alpha) .* (0.5 * D2 ./ K - alpha * log(K))
        else
            error("Unknown hyperparameter")
        end
    end

    return K
end
