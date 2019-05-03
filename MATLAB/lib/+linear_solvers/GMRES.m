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
        r = rhs - matrixMultiplyFunction(x);
        if(isPreconditioned)
            r = P \ r;
        end
        r_norm = norm(r);
        V(:,1) = r/r_norm;
        s(1) = r_norm;
        tempIter = 0;
        for k = 1:m
            num_iterations = num_iterations+1;

            % arnoldi
            [H, V] = arnold(H, V, k, n, matrixMultiplyFunction, isPreconditioned, P);

            % givens rotation
            [H, cs, sn, s] = apply_givens_rotation(H, cs, sn, s, i);
            resid = abs(s(k+1))/rhs_norm;
            if(resid < tolerance)
                tempIter = k;
                break;
            end
        end

        if (tempIter = 0)
            k = m;
        else
            k = tempIter;
        end
        y = H(1:k, 1:k) \ s(1:k);
        x = x + V(:, 1:k)*y;
        if(resid < tolerance)
            break;
        end
    end
end

function [H, V] = arnoldi(H, V, i, n, matrixMultiplyFunction, isPreconditioned, P)
    w = matrixMultiplyFunction(V(:,i));
    if(isPreconditioned)
        w = P\w;
    end
    for k = 1:i
        H(k, i) = w'*V(:,k);
        w = w - H(k, i)*V(:,k);
    end
    H(i+1,i, norm(w));
    V(:, i+1) = w/H(i+1,i);
end

function [H, cs, sn, s] apply_givens_rotation(H, cs, sn, s, i)
    for k = 1:i-1
        temp = cs(k)*H(k, i) + sn(k)*H(k+1, i);
        H(k+1, i) = -sn(k)*H(k, i) + cs(k)*H(k+1, i);
        H(k, i) = temp;
    end

    v1 = H(i, i);
    v2 = H(i+1, i);
    if(abs(v1) < eps)
        cs(i) = 0;
        sn(i) = 1;
    else
        t = sqrt(v1^2 + v2^2);
        cs(i) = abs(v1)/t;
        sn(i) = cs(i)*v2/v1;
    end

    s(i+1) = -sn(i)*s(i);
    s(i) = cs(i)*s(i);

    H(i, i) = cs(i)*H(i, i) + sn(i)*H(i+1, i);
    H(i+1, 1) = 0.0;
end
