function [ J ] = getFDThinFilmJacobian( q_FD, deltaX )
    num_cells = length(q_FD);
    odx4 = 1/deltaX^4;
    
    q3_iphalf = (0.5*(q_FD + q_FD([2:end,1]))).^3;
    q2_iphalf = (0.5*(q_FD + q_FD([2:end,1]))).^2;
    qxxx_iphalf = -q_FD([end,1:end-1]) + 3*q_FD - 3*q_FD([2:end,1]) + q_FD([3:end,1,2]);
    
    % This is the finite difference scheme
    % \dot(Q_i) = -(1/dx^4)*((1/2)(Q_i^3 + Q_{i+1}^3)*(-Q_{i-1} + 3*Q_i - 3*Q_{i+1} + Q_{i+2})
    % - ((1/2)(Q_{i-1}^3 + Q_{i}^3)(-Q_{i-2} + 3Q_{i-1} - 3Q_i + Q_{i+1})
    
    % a(i) = dG^i/dQ_{i-2}
    % a(i) = -1/dx^4 Q^3_{i-1/2}
    a = -odx4*q3_iphalf([end,1:end-1]);

    % b(i) = dG^i/dQ_{i-1}
    % b(i) = -1/dx^4*(-Q^3_{i+1/2} - 3*Q^3_{i-1/2} -
    % 3/2*Q^2_{i-1/2}*(Q_{xxx})_{i-1/2}
    b = -odx4*(-1.0*q3_iphalf - 3.0*q3_iphalf([end,1:end-1]) - 1.5*q2_iphalf([end,1:end-1]).*qxxx_iphalf([end,1:end-1]));
    
    % c(i) = dF^i/dQ_i
    % c(i) = -1/dx^4*(3Q^3_{i+1/2} + 3/2*Q^2_{i+1/2}*(Q_{xxx})_{i+1/2} + 3*Q^3_{i-1/2} - 3/2*Q^2_{i-1/2}*(Q_{xxx})_{i-1/2})
    c = -odx4*(3.0*q3_iphalf + 1.5*q2_iphalf.*qxxx_iphalf + 3.0*q3_iphalf([end,1:end-1]) - 1.5*q2_iphalf([end,1:end-1]).*qxxx_iphalf([end,1:end-1]));
    
    % d(i) = dF^i/dQ_{i+1}
    % d(i) = -1/dx^4*(-3*Q^3_{i+1/2} + 3/2*Q^2_{i+1/2}*(Q_{xxx})_{i+1/2} - Q^3_{i-1/2})
    d = -odx4*(-3.0*q3_iphalf + 1.5*q2_iphalf.*qxxx_iphalf - q3_iphalf([end,1:end-1]));
    
    % g(i) = dF^i/dQ_{i+2}
    % g(i) = -1/dx^4*f_{i+1/2} = -(1/2)(1/dx^4)(Q_i^3 + Q_{i+1}^3)
    g = -odx4*(q3_iphalf);
    
    % J(i,:) = [0, a(i),  b(i), c(i), d(i), g(i), 0]
    J = spdiags([[a(3:end);0;0], [b(2:end); 0], c, [0;d(1:end-1)], [0;0;g(1:end-2)]], -2:2, num_cells, num_cells);

    J(1,end-1) = a(1);
    J(1,end) = b(1);
    J(2,end) = a(2);
    J(end-1,1) = g(end-1);
    J(end,1) = d(end);
    J(end,2) = g(end);
end

