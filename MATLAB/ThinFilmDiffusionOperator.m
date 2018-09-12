function [Aq] = ThinFilmDiffusionOperator(q, dx, original_q, diffusivity, bc)
    [num_elems, num_basis_cpts] = size(q);

    odx = 1.0/dx;

    FQ = zeros(num_elems+4, 1);
    FR = zeros(num_elems+3, 1);
    FSRight = zeros(num_elems+2, 1);
    FSLeft  = zeros(num_elems+1, 1);
    etaLeft = zeros(num_elems+2, 1);
    etaRight = zeros(num_elems+2, 1);
    FU = zeros(num_elems+1, 1);

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
            % cell M - 1 = index num_elems+1
            FQ(1) = FQ(num_elems+1);
            FQ(2) = FQ(num_elems+2);
            FQ(end) = FQ(4);
            FQ(end-1) = FQ(3);
        case 'extrapolation'
            FQ(1) = FQ(3);
            FQ(2) = FQ(3);
            FQ(end-1) = FQ(num_elems+2);
            FQ(end) = FQ(num_elems+2);
    end

    % gradient, first derivative = R
    % R = 2/dx * q_{\xi}
    % R_i^k = 1/dx*(-int_{-1}{1}{Q \phi^k_{xi}(xi)}{xi} + FQ_{i+1/2}*phi^k(1) - FQ_{i-1/2}*phi^k(-1))
    switch bc
        case 'periodic'
            qExtended = [q(end-1,:); q(end,:); q; q(1,:)];
        case 'extrapolation'
            qExtended = [q(1,:); q(1,:); q; q(end,:)];
    end
    R = odx*(-(Smat*qExtended')' + FQ(2:end)*phi - FQ(1:end-1)*(sgn.*phi));

    % fluxes to evaluate second derivative
    % FR always evaluated on left side
    % FR[i] = Sum_{l = 1}^{num_basis_cpts}{R_{i-1}^l *phi^l(1)}
    % FR[i] = R evaluated at (i - 1/2) boundary
    FR(:) = phi*R';

    % Second derivative - S
    % S = 2/dx r_{\xi} = 4/dx^2 q_{\xi, \xi}
    S = odx*(-(Smat*R(2:end,:)')' + FR(2:end)*phi - FR(1:end-1)*(sgn.*phi));

    % fluxes to evaluate U
    % SRight[i] = Sum_{l = 1}^{num_basis_cpts}{S_{i}^l *phi^l(-1)}
    % SRight[i] = S evaluated on right side of (i - 1/2) boundary
    % SLeft[i] = Sum_{l = 1}^{num_basis_cpts}{S_{i}^l *phi^l(1)}
    % SLeft[i] = S evaluated on left side of (i + 1/2) boundary
    % etaLeft[i] = eta evaluated on left side of (i-1/2) boundary
    % etaLeft[i] = eta_{i}(1)
    % etaRight[i] = eta evaluate on right side of (i-1/2) boundary
    % etaRight[i] = eta_i(-1)
    FSRight(:) = (sgn.*phi)*S';
    FSLeft(:)    = phi*S(1:end-1)';
    etaRight(2:end-1) = ((sgn.*phi)*original_q').^3;
    etaLeft(3:end)    = (phi*original_q').^3;
    switch bc
        case 'periodic'
            etaRight(1) = etaRight(end-1);
            etaRight(end) = etaRight(2);
            etaLeft(1) = etaLeft(end-1);
            etaLeft(2) = etaLeft(end);
        case 'extrapolation'
            etaRight(1) = etaRight(2);
            etaRight(end) = etaRight(end-1);
            etaLeft(1) = etaLeft(3);
            etaLeft(2) = etaLeft(3);
    end
    etaAverage = (etaRight + etaLeft)/2;

    % Q^3 * Third Derivate - U
    Stmp = zeros(num_elems+1,num_basis_cpts);
    eta = @(i, xi) sum(arrayfun(@(l) original_q(i,l)*legendrePolynomial(l,xi), 1:num_basis_cpts)).^3;
    for i = 1:num_elems+1
        % bc for eta
        im1 = i-1;
        if(i==1)
            switch bc
                case 'periodic'
                    im1 = num_elems;
                case 'extrapolation'
                    im1 = 1;
            end
        end
        for k = 1:num_basis_cpts
            integrandFunction = @(xi, l) legendrePolynomial(k, xi)*legendrePolynomialDerivative(l, xi)*eta(im1, xi);
            Stmp(i,k) = sum(arrayfun(@(l) S(i, l)*gaussQuadrature(num_basis_cpts, @(xi) integrandFunction(xi, l)), 2:num_basis_cpts));
        end
    end
    U = odx*(Stmp + etaAverage(2:end)*phi.*SRight(2:end) - etaLeft(2:end)*phi.*SLeft(1:end) + ...
        (etaRight(1:end-1) + etaAverage(1:end-1))*(sgn.*phi).*SRight(1:end-1));

    % fluxes to evaluate solution
    % FU[i] is flux of U at i - 1/2 interface = U_{i-1}(x_{i-1/2}) = flux evaluated on left hand side
    % FU[i] = Sum_{l = 1}^{num_basis_cpts}{U_{i-1}^l *phi^l(1)}
    FU(:) = phi*U';

    Aq = diffusivity*odx*((Smat*U')' - FU(2:end)*phi + FU(1:end-1)*(sgn.*phi));
end
