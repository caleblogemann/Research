function [A] = GetDiscretizationMatrix(q_old, q_exact_function, t, dt, dx, forcing_function)
    % approximate q(t + dt/2) as q(t) + dt/2*q_t(t) = q(t) + dt/2 *((-q^3q_xxx)_x|t + f(t))
    q_exact_t = q_exact_function(t);
    A = getFDThinFilmMatrix(q_exact_t, dx);
    q_exact_t_FD = dog_math.getQFD(q_exact_t);
    mq3qxxx_x = A*q_exact_t_FD;

    q_t_half_dt = q_exact_t;
    q_t_half_dt(:, 1, 1) = q_t_half_dt(:, 1, 1) + 0.5*dt*(mq3qxxx_x + forcing_function(t));
    A = getFDThinFilmMatrix(q_t_half_dt, dx);
end
