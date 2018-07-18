nBasisCpts = 2;
nCells = 10;
a = 0;
b = 1;
qleft = 0.3323;
qright = 0.1;
discontinuity = 0.5;
deltaX = (b - a)/nCells;
riemannFunction = @(x) qleft*double(x < discontinuity) + qright*double(x >= discontinuity);
q = projectQ(riemannFunction, nBasisCpts, nCells, a, b);
[LDG] = getLDGThinFilmMatrix(q, nCells, nBasisCpts, deltaX);
bc = 'extrapolation';
[num_elems, num_basis_cpts] = size(q);

    odx = 1.0/dx;

    R = zeros(num_elems,num_basis_cpts);
    S = zeros(num_elems, num_basis_cpts);
    U = zeros(num_elems, num_basis_cpts);

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