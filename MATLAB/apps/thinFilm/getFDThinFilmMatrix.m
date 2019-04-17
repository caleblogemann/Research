function [A] = getFDThinFilmMatrix(q_FD, deltaX)
    % q_DG is a dg solution used to generate F = q^3, q_DG is the q before cubing
    % assume for now periodic boundary conditions
    % also use diffusivity constant 1

    %[num_cells, num_eqns, num_basis_cpts] = size(q_DG);
    % num_eqns should equal 1, remove to be able to multiply with q_DG
    %q_DG = squeeze(q_DG);
    num_cells = length(q_FD);
    odx4 = 1/deltaX^4;

    %phi1 = basis.getPhiVector(1, num_basis_cpts);
    %phim1 = basis.getPhiVector(-1, num_basis_cpts);

    % f(i) = 0.5*((Q_{i-1/2}^-)^3 + (Q_{i-1/2}^+)^3)
    % f(1) = f_{1/2}
    % f(i) = f_{i-1/2}
    % f(end) = f(nCells) = f_{N-1/2}
    %f = (0.5*((q_DG([end, 1:end-1], :)*phi1) + (q_DG*phim1))).^3;
    f = (0.5*(q_FD([end,1:end-1]) + q_FD)).^3;
    fimhalf = f;
    fiphalf = f([2:end, 1]);

    % FD = a*Q_{i-2} + b*Q_{i-1} + c*Q_i + d*Q_{i+1} + g*Q_{i+2}
    a = -1.0*fimhalf*odx4;
    b = (fiphalf + 3.0*fimhalf)*odx4;
    c = -3.0*(fiphalf + fimhalf)*odx4;
    d = (3.0*fiphalf + fimhalf)*odx4;
    g = -1.0*fiphalf*odx4;

    A = spdiags([[a(3:end);0;0], [b(2:end); 0], c, [0;d(1:end-1)], [0;0;g(1:end-2)]], -2:2, num_cells, num_cells);

    % periodic boundary conditions
    A(1,end-1) = a(1);
    A(1,end) = b(1);
    A(2,end) = a(2);
    A(end-1,1) = g(end-1);
    A(end,1) = d(end);
    A(end,2) = g(end);
end
