function [ G ] = FDThinFilmOperator( q_FD, deltaX )
  	num_cells = length(q_FD);
    odx4 = 1/deltaX^4;

    % q3_iphalf(i) = (0.5*(Q_{i} + Q_{i+1}))^3
    q3_iphalf = (0.5*(q_FD + q_FD([2:end,1]))).^3;

    qxxx_iphalf = (-q_FD([end,1:end-1]) + 3*q_FD - 3*q_FD([2:end,1]) + q_FD([3:end,1,2]));

    % G = -1/dx4*(Q^3_{i+1/2}*(Q_{xxx})_{i+1/2} - Q^3_{i-1/2}*(Q_{xxx})_{i-1/2})
    G = -odx4*(q3_iphalf.*qxxx_iphalf - q3_iphalf([end,1:end-1]).*qxxx_iphalf([end,1:end-1])); 
end

