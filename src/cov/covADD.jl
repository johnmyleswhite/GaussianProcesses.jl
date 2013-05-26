function covADD(cov::Vector, hyp::Vector, x::Matrix, z::Matrix, i)
    # Additive covariance function using a 1d base covariance function 
    # cov(x^p, x^q; hyp) with individual hyperparameters hyp.
    #
    # k(x^p, x^q) = \sum_{r \in R} sf_r \sum_{|I|=r}
    #                 \prod_{i \in I} cov(x^p_i, x^q_i; hyp_i)
    #
    # hyp = [hyp_1
    #        hyp_2
    #         ...
    #        hyp_D 
    #        log(sf_R(1))
    #         ...
    #        log(sf_R(end))]
    #
    # where hyp_d are the parameters of the 1d covariance function which are shared
    # over the different values of R(1) to R(end).
    #
    # Please see the paper Additive Gaussian Processes by Duvenaud, Nickisch and 
    # Rasmussen, NIPS, 2011 for details.
    #
    # Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
    #
    # See also COVFUNCTIONS.M.

    R = cov[1]

    # number of hypers per individual covariance
    nh = eval(feval(cov[2]))

    # number of different degrees of interaction
    nr = length(R)

    # report number of hyper parameters
    if nargin < 3
        K = ["D*", int2str(nh), "+", int2str(nr)]
        return K
    end

    # make sure, z exists
    if nargin < 4
        z = []
    end

    # determine mode
    xeqz = length(z) == 0
    dg = strcmp(z, "diag") && length(z) > 0

    # dimensionality
    n, D = size(x)

    # signal variances of individual degrees
    sf2 = exp(2 * hyp[D * nh + 1:nr])

    # evaluate dimensionwise covariances K
    Kd = Kdim(cov[2], hyp, x, z)

    # covariances
    if nargin < 5
        # Rth elementary symmetric polynomials
        EE = elsympol(Kd, max(R))
        K = 0
        # sf2 weighted sum
        for ii in 1:nr
            K = K + sf2[ii] * EE[:, :, R[ii] + 1]
        end
    else
        # derivatives
        if i <= D * nh
            # individual covariance function parameters
            # j is the dimension of the hyperparameter
            j = fix(1 + (i - 1) / nh)
            if dg
                zj = "diag"
            else
                if xeqz
                    zj = []
                else
                    zj = z[:, j]
                end
            end

            # other dK = 0
            dKj = feval(cov[2],
                        hyp[nh * (j - 1) + 1:nh],
                        x[:, j],
                        zj,
                        i - (j - 1) * nh)

            # the final derivative is a sum of multilinear terms,
            # so if only one term depends on the hyperparameter under
            # consideration, we can factorise it  out and compute the sum
            # with one degree less

            # R-1th elementary sym polyn
            E = elsympol(Kd[:, :, [1:(j - 1), (j + 1):D]], max(R) - 1)

            # sf2 weighted sum
            K = 0
            for ii in 1:nr
                K = K + sf2[ii] * E[:, :, R[ii]]
            end
            K = dKj .* K
        elseif i <= D * nh + nr
            # Rth elementary symmetric polynomials
            EE = elsympol(Kd, max(R))
            j = i - D * nh

            # rest of the sf2 weighted sum
            K = 2 * sf2[j] * EE[:, :, R[j] + 1]
        else
            error("Unknown hyperparameter")
        end
    end
end

# evaluate dimensionwise covariances K
function Kdim(cov::Vector, hyp::Vector, x::Matrix, z::Matrix)
    # dimensionality
    n, D = size(x)

    # number of hypers per individual covariance
    nh = eval(feval(cov))

    # make sure, z exists
    if nargin < 4
        z = []
    end

    # determine mode
    xeqz = length(z) == 0
    dg = strcmp(z, "diag") && length(z) > 0

    # allocate memory
    if dg
        K = zeros(n, 1, D)
    else
        if xeqz
            K = zeros(n, n, D)
        else
            K = zeros(n, size(z, 1), D)
        end
    end

    for d in 1:D                               
        # hyperparameter of dimension d
        hyp_d = hyp[nh * (d - 1) + 1:nh]
        if dg
            K[:, :, d] = feval(cov, hyp_d, x[:, d], "diag")
        else
            if xeqz
                K[:, :, d] = feval(cov, hyp_d, x[:, d])
            else
                K[:, :, d] = feval(cov, hyp_d, x[:, d], z[:, d])
            end
        end
    end

    return K
end
