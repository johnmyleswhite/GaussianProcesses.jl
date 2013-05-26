function covFITC(cov::Vector, xu, hyp::Vector, x::Matrix, z::Matrix, i)
    # Covariance function to be used together with the FITC approximation
    #
    # The function allows for more than one output argument and does not respect the
    # interface of a proper covariance function. In fact, it wraps a proper
    # covariance function such that it can be used together with infFITC.m.
    # Instead of outputing the full covariance, it returns cross-covariances between
    # the inputs x, z and the inducing inputs xu as needed by infFITC.m
    #
    # Copyright (c) by Ed Snelson, Carl Edward Rasmussen and Hannes Nickisch, 2011-11-02.
    #
    # See also COVFUNCTIONS.M, INFFITC.M.

    if nargin < 4
        K = feval(cov[:])
        return
    end

    # make sure, z exists
    if nargin < 5
        z = []
    end

    # determine mode
    xeqz = length(z) == 0
    dg = strcmp(z, "diag") && length(z) > 0

    if size(xu, 2) != size(x, 2)
        error("Dimensionality of inducing inputs must match training inputs")
    end

    # covariances
    if nargin < 6
        if dg
            K = feval(cov[:], hyp, x, "diag")
        else
            if xeqz
                K = feval(cov[:], hyp, x, "diag")
                if nargout > 1
                    Kuu = feval(cov[:], hyp, xu)
                end
                if nargout > 2
                    Ku = feval(cov[:], hyp, xu, x)
                end
            else
                K = feval(cov[:], hyp, xu, z)
            end
        end
    # derivatives
    else
        if dg
            K = feval(cov[:], hyp, x, "diag", i)
        else
            if xeqz
                K = feval(cov[:], hyp, x, "diag", i)
                if nargout > 1
                    Kuu = feval(cov[:], hyp, xu, [], i)
                end
                if nargout > 2
                    Ku = feval(cov[:], hyp, xu, x, i)
                end
            else
                K = feval(cov[:], hyp, xu, z, i)
            end
        end
    end

    return K, Kuu, Ku
end
