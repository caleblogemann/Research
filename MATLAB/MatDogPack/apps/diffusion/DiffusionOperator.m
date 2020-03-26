function [Aq] = DiffusionOperator(q, dx, diffusivity, bc)
    [num_elems, num_basis_cpts] = size(q);

    odx = 1.0/dx;

    FQ = zeros(num_elems+2, 1);
    FR = zeros(num_elems+1, 1);

    phi = arrayfun(@(k) legendrePolynomial(k, 1), 1:num_basis_cpts);

    sgn = [1.0, -1.0, 1.0, -1.0, 1.0];
    sgn = sgn(1:num_basis_cpts);
    % Smat[k][l] represents int_{-1}^{1} \phi^{k}_{\xi}{\xi} \phi^{l}(\xi) d\xi
    Smat = [0.0, 0.0, 0.0, 0.0, 0.0;
        2.0*sqrt(3), 0.0, 0.0, 0.0, 0.0;
        0.0, 2.0*sqrt(3)*sqrt(5), 0.0, 0.0, 0.0;
        2.0*sqrt(7), 0.0, 2.0*sqrt(5)*sqrt(7), 0.0, 0.0;
        0.0, 6.0*sqrt(3), 0.0, 6.0*sqrt(7), 0.0];
    Smat = Smat(1:num_basis_cpts, 1:num_basis_cpts);

    % fluxes to evaluate gradient
    % Q flux always evaluated at right side of boundary
    % FQ[i] = Sum_{l = 1}^{num_basis_cpts}{Q_i^l *phi^l(-1)}
    % FQ[i] = Q at boundary (i - 3/2)
    % FQ[1] = Q at boundary (-1/2)
    FQ(2:num_elems+1) = (sgn.*phi)*q';
    switch bc
        case 'periodic'
            FQ(1) = FQ(end-1);
            FQ(end) = FQ(2);
        case 'extrapolation'
            FQ(1) = FQ(2);
            FQ(end) = FQ(end-1);
    end

    % gradient, first derivative = R
    % R = 2/dx * q_{\xi}
    % R_i^k = 1/dx*(-int_{-1}{1}{Q \phi^k_{xi}(xi)}{xi} + FQ_{i+1/2}*phi^k(1) - FQ_{i-1/2}*phi^k(-1))
    switch bc
        case 'periodic'
            qExtended = [q(end,:); q];
        case 'extrapolation'
            qExtended = [q(1,:); q];
    end
    R = odx*(-(Smat*qExtended')' + FQ(2:end)*phi - FQ(1:end-1)*(sgn.*phi));

    % fluxes to evaluate second derivative
    % FR always evaluated on left side
    % FR[i] = Sum_{l = 1}^{num_basis_cpts}{R_{i-1}^l *phi^l(1)}
    % FR[i] = R evaluated at (i - 1/2) boundary
    FR(:) = phi*R';

    % Second derivative - Aq, Q
    % Aq = 2/dx r_{\xi} = 4/dx^2 q_{\xi, \xi}
    Aq = diffusivity*odx*(-(Smat*R(2:end, :)')' + FR(2:end)*phi - FR(1:end-1)*(sgn.*phi));
end
