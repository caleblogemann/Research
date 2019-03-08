function [x, resid, num_iterations] = GMRES(rhs, x, matrixMultiplyFunction, max_iterations, tolerance, hasRestarts, restartIterations, isPreconditioned, P)
    num_iterations = 0;

    rhs_norm = norm(rhs);
    if(rhs_norm == 0.0)
        rhs_norm = 1.0;
    end

    n = length(x);

    if(hasRestarts)
        m = restartIterations;
        max_number_restarts = ceil(max_iterations/m);
    else
        m = max_iterations;
        max_number_restarts = 1;
    end

    V = zeros(n, m+1);
    H = zeros(m+1, m);
    cs = zeros(m, 1);
    sn = zeros(m, 1);

    s = zeros(m+1, 1);

    y = zeros(m, 1);

    Vy = zeros(n, 1);

    for restartIter = 1:max_number_restarts
        r = rhs - A*x;
        if(isPreconditioned)
            r = P \ r;
        end
        r_norm = norm(r);
        V(:,1) = r/r_norm;
        s(1) = r_norm;
        for k = 1:m
            num_iterations = num_iterations+1;
        end

    end

end
