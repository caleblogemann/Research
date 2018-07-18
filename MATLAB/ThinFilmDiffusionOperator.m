function [Aq] = ThinFilmDiffusionOperator(q, dx, original_q, bc)
    [num_elems, num_basis_cpts] = size(q);

    odx = 1.0/dx;

    %R = zeros(num_elems,num_basis_cpts);
    %S = zeros(num_elems, num_basis_cpts);
    %U = zeros(num_elems, num_basis_cpts);

    FQ = zeros(num_elems+1, 1);
    FR = zeros(num_elems+1, 1);
    SRight = zeros(num_elems+1, 1);
    SLeft  = zeros(num_elems+1, 1);
    FU = zeros(num_elems+1, 1);
    etaLeft = zeros(num_elems+1, 1);
    etaRight = zeros(num_elems+1, 1);

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
    % FQ[i] = Q at boundary (i - 1/2)
    FQ(1:num_elems) = (sgn.*phi)*q';
    switch bc
        case 'periodic'
            FQ(end) = FQ(1);
        case 'extrapolation'
            FQ(end) = FQ(end-1);
    end

    % gradient, first derivative = R
    % R = 2/dx * q_{\xi}
    % R_i^k = 1/dx*(-int_{-1}{1}{Q \phi^k_{xi}(xi)}{xi} + FQ_{i+1/2}*phi^k(1) - FQ_{i-1/2}*phi^k(-1))
    R = odx*(-(Smat*q')' + FQ(2:end)*phi - FQ(1:end-1)*(sgn.*phi));

    % fluxes to evaluate second derivative
    % FR always evaluated on left side
    % FR[i] = Sum_{l = 1}^{num_basis_cpts}{R_{i-1}^l *phi^l(1)}
    % FR[i] = R evaluated at (i - 1/2) boundary
    FR(2:end) = phi*R';
    switch bc
        case 'periodic'
            FR(1) = FR(end);
        case 'extrapolation'
            FR(1) = FR(2);
    end

    % Second derivative - S
    % S = 2/dx r_{\xi} = 4/dx^2 q_{\xi, \xi}
    S = odx*(-(Smat*R')' + FR(2:end)*phi - FR(1:end-1)*(sgn.*phi));

    % fluxes to evaluate U
    % SRight[i] = Sum_{l = 1}^{num_basis_cpts}{S_{i}^l *phi^l(-1)}
    % SRight[i] = S evaluated on right side of (i - 1/2) boundary
    % SLeft[i] = Sum_{l = 1}^{num_basis_cpts}{S_{i}^l *phi^l(1)}
    % SLeft[i] = S evaluated on left side of (i - 1/2) boundary
    % etaLeft[i] = eta evaluated on left side of (i-1/2) boundary
    % etaLeft[i] = eta_{i}(1)
    % etaRight[i] = eta evaluate on right side of (i-1/2) boundary
    % etaRight[i] = eta_i(-1)
    SRight(1:end-1) = (sgn.*phi)*S';
    SLeft(2:end)    = phi*S';
    etaRight(1:end-1) = ((sgn.*phi)*original_q').^3;
    etaLeft(2:end)    = (phi*original_q').^3;
    switch bc
        case 'periodic'
            SRight(end) = SRight(1);
            SLeft(1) = SLeft(end);
            etaRight(end) = etaRight(1);
            etaLeft(1) = etaLeft(end);
        case 'extrapolation'
            SRight(end) = SRight(end-1);
            SLeft(1) = SLeft(2);
            etaRight(end) = etaRight(end-1);
            etaLeft(1) = etaLeft(2);
    end

    % Q^3 * Third Derivate - U
    Stmp = zeros(num_elems,num_basis_cpts);
    eta = @(i, xi) sum(arrayfun(@(l) original_q(i,l)*legendrePolynomial(l,xi), 1:num_basis_cpts)).^3;
    for i = 1:num_elems
        for k = 1:num_basis_cpts
            integrandFunction = @(xi, l) legendrePolynomial(k, xi)*legendrePolynomialDerivative(l, xi)*eta(i, xi);
            Stmp(i,k) = sum(arrayfun(@(l) S(i, l)*gaussQuadrature(num_basis_cpts, @(xi) integrandFunction(xi, l)), 2:num_basis_cpts));
        end
    end
    U = odx*(Stmp + 0.5*(etaRight(2:end) + etaLeft(2:end))*phi.*SRight(2:end) - ...
        (etaLeft(2:end).*SLeft(2:end))*phi - (0.5*(etaLeft(1:end-1) - ...
        etaRight(1:end-1)).*SRight(1:end-1))*(sgn.*phi));

    % fluxes to evaluate solution
    % FU[i] is flux of U at i - 1/2 interface = U_{i-1}(x_{i-1/2}) = flux evaluated on left hand side
    % FU[i] = Sum_{l = 1}^{num_basis_cpts}{U_{i-1}^l *phi^l(1)}
    FU(2:end) = phi*U';
    switch bc
        case 'periodic'
            FU(1) = FU(end);
        case 'extrapolation'
            FU(1) = FU(2);
    end

    Aq = odx*((Smat*U')' - FU(2:end)*phi + FU(1:end-1)*(sgn.*phi));
end
