function covSEiso(hyp::Vector, x::Matrix, z::Matrix, i)
    # Squared Exponential covariance function with isotropic distance measure. The 
    # covariance function is parameterized as:
    #
    # k(x^p,x^q) = sf^2 * exp(-(x^p - x^q)'*inv(P)*(x^p - x^q)/2) 
    #
    # where the P matrix is ell^2 times the unit matrix and sf^2 is the signal
    # variance. The hyperparameters are:
    #
    # hyp = [ log(ell)
    #         log(sf)  ]
    #
    # For more help on design of covariance functions, try "help covFunctions".
    #
    # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
    #
    # See also COVFUNCTIONS.M.#

    # report number of parameters
    if nargin < 2
        K = "2"
        return K
    end

    # make sure, z exists
    if nargin < 3
        z = []
    end

    # determine mode
    xeqz = length(z) == 0
    dg = strcmp(z, "diag") && length(z) > 0

    # characteristic length scale
    ell = exp(hyp[1])
    # signal variance
    sf2 = exp(2 * hyp[2])

    # precompute squared distances
    if dg                                                               
        # vector kxx
        K = zeros(size(x, 1), 1)
    else
        if xeqz
            # symmetric matrix Kxx
            K = sq_dist(x' / ell)
        else
            # cross covariances Kxz
            K = sq_dist(x' / ell, z' / ell)
        end
    end

    # covariances
    if nargin < 4
        K = sf2 * exp(-K / 2)
    # derivatives
    else                                                               
        if i == 1
            K = sf2 * exp(-K / 2) .* K
        elseif i == 2
            K = 2 * sf2 * exp(-K / 2)
        else
            error("Unknown hyperparameter")
        end
    end

    return K
end
