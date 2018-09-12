a = 0;
b = 1;
qleft = 0.3323;
qright = 0.1;
discontinuity = 0.5;
riemannFunction = @(x) qleft*double(x < discontinuity) + qright*double(x >= discontinuity);
constantFunction = @(x) 1;
mu = 0.5;
sigma = 0.05;
height = 0.3;
gaussianFunction = @(x) height*exp(-1/2*((x - mu)/sigma).^2);
num_elems = 4;
deltaX = (b - a)/num_elems;
tolerance = 1e-8;
diffusivity = 1.0;

%% Diffusion LDG
    disp('Diffusion LDG');
    % Should be zero for constant function
        disp('Constant Function');
        % extrapolation
            bc = 'extrapolation';
            disp(bc);
            % 1 basis cpts
                disp('nBasisCpts = 1');
                num_basis_cpts = 1;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                disp(norm(A*q) <= tolerance);
            % 2 basis cpts
                disp('nBasisCpts = 2');
                num_basis_cpts = 2;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                disp(norm(A*getVector(q)) <= tolerance);
        % periodic
            bc = 'periodic';
            disp(bc);
            % 1 basis cpts
                disp('nBasisCpts = 1');
                num_basis_cpts = 1;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                disp(norm(A*q) <= tolerance)
            % 2 basis cpts
                disp('nBasisCpts = 2');
                num_basis_cpts = 2;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                disp(norm(A*getVector(q)) <= tolerance)

    % Should match FD version for extrapolation num_basis_cpts = 1, diffusivity = 1.0,
        disp('Matching FD Matrix');
        num_basis_cpts = 1;
        bc = 'extrapolation';
        disp(bc);
        ALDG = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
        AFD = getFDDiffusionMatrix(num_elems, diffusivity, deltaX, bc);
        disp(norm(ALDG - AFD) <= tolerance)

    % and for periodic
        bc = 'periodic';
        disp(bc);
        ALDG = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
        AFD = getFDDiffusionMatrix(num_elems, diffusivity, deltaX, bc);
        disp(norm(ALDG - AFD) <= tolerance)

    % Should match operator version with extrapolation
        disp('Matching Operator version');
        bc = 'extrapolation';
        disp(bc);
        % at num_basis_cpts = 1
            num_basis_cpts = 1;
            disp('nBasisCpts = 1');
            ALDG = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
            Aq = DiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - Aq) <= tolerance)

        % at num_basis_cpts = 2
            disp('nBasisCpts = 2');
            num_basis_cpts = 2;
            ALDG = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
            Aq = DiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - getVector(Aq)) <= tolerance)

        % at num_basis_cpts = 3
            disp('nBasisCpts = 3');
            num_basis_cpts = 3;
            ALDG = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
            Aq = DiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - getVector(Aq)) <= tolerance)

    % And with periodic boundary conditions
        bc = 'periodic';
        disp(bc);
        % at num_basis_cpts = 1
            disp('nBasisCpts = 1');
            num_basis_cpts = 1;
            ALDG = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(gaussianFunction, num_basis_cpts, num_elems, a, b);
            Aq = DiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - getVector(Aq)) <= tolerance)

        % at num_basis_cpts = 2
            disp('nBasisCpts = 2');
            num_basis_cpts = 2;
            ALDG = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(gaussianFunction, num_basis_cpts, num_elems, a, b);
            Aq = DiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - getVector(Aq)) <= tolerance)

        % at num_basis_cpts = 3
            disp('nBasisCpts = 3');
            ALDG = getLDGDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(gaussianFunction, num_basis_cpts, num_elems, a, b);
            Aq = DiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - getVector(Aq)) <= tolerance)

%% HyperDiffusion
    disp('Hyper Diffusion LDG');
    num_elems = 6;
    deltaX = (b - a)/num_elems;
    % Should be zero for constant function
        disp('Constant Function');
        % extrapolation
            bc = 'extrapolation';
            disp(bc);
            % 1 basis cpts
                disp('nBasisCpts = 1');
                num_basis_cpts = 1;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                disp(norm(A*q) <= tolerance);
            % 2 basis cpts
                disp('nBasisCpts = 2');
                num_basis_cpts = 2;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                disp(norm(A*getVector(q)) <= tolerance);
        % periodic
            bc = 'periodic';
            disp(bc);
            % 1 basis cpts
                disp('nBasisCpts = 1');
                num_basis_cpts = 1;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                disp(norm(A*q) <= tolerance)
            % 2 basis cpts
                disp('nBasisCpts = 2');
                num_basis_cpts = 2;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                disp(norm(A*getVector(q)) <= tolerance)

    % Should match FD version for extrapolation num_basis_cpts = 1, diffusivity = 1.0,
        disp('Matching FD Matrix');
        num_basis_cpts = 1;
        bc = 'extrapolation';
        disp(bc);
        ALDG = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
        AFD = getFDHyperDiffusionMatrix(num_elems, diffusivity, deltaX, bc);
        disp(norm(ALDG - AFD) <= tolerance)

    % and for periodic
        bc = 'periodic';
        disp(bc);
        ALDG = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
        AFD = getFDHyperDiffusionMatrix(num_elems, diffusivity, deltaX, bc);
        disp(norm(ALDG - AFD) <= tolerance)

    % Should match operator version with extrapolation
        disp('Matching Operator version');
        bc = 'extrapolation';
        disp(bc);
        % at num_basis_cpts = 1
            num_basis_cpts = 1;
            disp('nBasisCpts = 1');
            ALDG = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
            Aq = HyperDiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - Aq) <= tolerance)

        % at num_basis_cpts = 2
            disp('nBasisCpts = 2');
            num_basis_cpts = 2;
            ALDG = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
            Aq = HyperDiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - getVector(Aq)) <= tolerance)

        % at num_basis_cpts = 3
            disp('nBasisCpts = 3');
            num_basis_cpts = 3;
            ALDG = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
            Aq = HyperDiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - getVector(Aq)) <= tolerance)

    % And with periodic boundary conditions
        bc = 'periodic';
        disp(bc);
        % at num_basis_cpts = 1
            disp('nBasisCpts = 1');
            num_basis_cpts = 1;
            ALDG = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(gaussianFunction, num_basis_cpts, num_elems, a, b);
            Aq = HyperDiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - getVector(Aq)) <= tolerance)

        % at num_basis_cpts = 2
            disp('nBasisCpts = 2');
            num_basis_cpts = 2;
            ALDG = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(gaussianFunction, num_basis_cpts, num_elems, a, b);
            Aq = HyperDiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - getVector(Aq)) <= tolerance)

        % at num_basis_cpts = 3
            disp('nBasisCpts = 3');
            ALDG = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
            q = projectQ(gaussianFunction, num_basis_cpts, num_elems, a, b);
            Aq = HyperDiffusionOperator(q, deltaX, diffusivity, bc);
            qVector = getVector(q);
            disp(norm(ALDG*qVector - getVector(Aq)) <= tolerance)

%% Thin Film Diffusion
    disp('Thin Film Diffusion LDG');
    num_elems = 6;
    deltaX = (b - a)/num_elems;
    % Should be zero for constant function
        disp('Constant Function');
        % extrapolation
            bc = 'extrapolation';
            disp(bc);
            % 1 basis cpts
                disp('nBasisCpts = 1');
                num_basis_cpts = 1;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGThinFilmMatrix(q, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, q, diffusivity, bc);
                disp(norm(A*q) <= tolerance);
                disp(norm(Aq) <= tolerance);
            % 2 basis cpts
                disp('nBasisCpts = 2');
                num_basis_cpts = 2;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGThinFilmMatrix(q, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, q, diffusivity, bc);
                disp(norm(A*getVector(q)) <= tolerance);
                disp(norm(Aq) <= tolerance);

        % periodic
            bc = 'periodic';
            disp(bc);
            % 1 basis cpts
                disp('nBasisCpts = 1');
                num_basis_cpts = 1;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGThinFilmMatrix(q, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, q, diffusivity, bc);
                disp(norm(A*q) <= tolerance);
                disp(norm(Aq) <= tolerance);
            % 2 basis cpts
                disp('nBasisCpts = 2');
                num_basis_cpts = 2;
                q = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGThinFilmMatrix(q, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, q, diffusivity, bc);
                disp(norm(A*getVector(q)) <= tolerance);
                disp(norm(Aq) <= tolerance);

    % Should Match hyper diffusion when original_q is one
        disp('Matching HyperDiffusion');
        % extrapolation
            bc = 'extrapolation';
            disp(bc);
            % 1 basis cpts
                disp('nBasisCpts = 1');
                num_basis_cpts = 1;
                one = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                q = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGThinFilmMatrix(one, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Hyper = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, one, diffusivity, bc);
                disp(norm(A - Hyper) <= tolerance);
                disp(norm(Aq - Hyper*q) <= tolerance);
            % 2 basis cpts
                disp('nBasisCpts = 2');
                num_basis_cpts = 2;
                one = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                q = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGThinFilmMatrix(one, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Hyper = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, one, diffusivity, bc);
                disp(norm(A - Hyper) <= tolerance);
                disp(norm(getVector(Aq) - Hyper*getVector(q)) <= tolerance);

        % periodic
            bc = 'periodic';
            disp(bc);
            % 1 basis cpts
                disp('nBasisCpts = 1');
                num_basis_cpts = 1;
                one = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                q = projectQ(gaussianFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGThinFilmMatrix(one, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Hyper = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, one, diffusivity, bc);
                disp(norm(A - Hyper) <= tolerance);
                disp(norm(Aq - Hyper*q) <= tolerance);
            % 2 basis cpts
                disp('nBasisCpts = 2');
                num_basis_cpts = 2;
                one = projectQ(constantFunction, num_basis_cpts, num_elems, a, b);
                q = projectQ(gaussianFunction, num_basis_cpts, num_elems, a, b);
                A = getLDGThinFilmMatrix(one, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Hyper = getLDGHyperDiffusionMatrix(num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, one, diffusivity, bc);
                disp(norm(A - Hyper) <= tolerance);
                disp(norm(getVector(Aq) - Hyper*getVector(q)) <= tolerance);

    % Should match operator version
        disp('Matching Operator version');
        % with extrapolation
            bc = 'extrapolation';
            disp(bc);
            % at num_basis_cpts = 1
                num_basis_cpts = 1;
                disp('nBasisCpts = 1');
                q = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
                ALDG = getLDGThinFilmMatrix(q, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, q, diffusivity, bc);
                qVector = getVector(q);
                disp(norm(ALDG*qVector - Aq) <= tolerance)
            % at num_basis_cpts = 2
                disp('nBasisCpts = 2');
                num_basis_cpts = 2;
                q = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
                ALDG = getLDGThinFilmMatrix(q, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, q, diffusivity, bc);
                qVector = getVector(q);
                disp(norm(ALDG*qVector - Aq) <= tolerance)
            % at num_basis_cpts = 3
                disp('nBasisCpts = 3');
                num_basis_cpts = 3;
                q = projectQ(riemannFunction, num_basis_cpts, num_elems, a, b);
                ALDG = getLDGThinFilmMatrix(q, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, q, diffusivity, bc);
                qVector = getVector(q);
                disp(norm(ALDG*qVector - Aq) <= tolerance)

        % And with periodic boundary conditions
            bc = 'periodic';
            disp(bc);
            % at num_basis_cpts = 1
                disp('nBasisCpts = 1');
                num_basis_cpts = 1;
                q = projectQ(gaussianFunction, num_basis_cpts, num_elems, a, b);
                ALDG = getLDGThinFilmMatrix(q, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, q, diffusivity, bc);
                qVector = getVector(q);
                disp(norm(ALDG*qVector - Aq) <= tolerance)
            % at num_basis_cpts = 2
                disp('nBasisCpts = 2');
                num_basis_cpts = 2;
                q = projectQ(gaussianFunction, num_basis_cpts, num_elems, a, b);
                ALDG = getLDGThinFilmMatrix(q, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, q, diffusivity, bc);
                qVector = getVector(q);
                disp(norm(ALDG*qVector - Aq) <= tolerance)
            % at num_basis_cpts = 3
                disp('nBasisCpts = 3');
                num_basis_cpts = 3;
                q = projectQ(gaussianFunction, num_basis_cpts, num_elems, a, b);
                ALDG = getLDGThinFilmMatrix(q, num_elems, num_basis_cpts, deltaX, diffusivity, bc);
                Aq = ThinFilmDiffusionOperator(q, deltaX, q, diffusivity, bc);
                qVector = getVector(q);
                disp(norm(ALDG*qVector - Aq) <= tolerance)

%% Thin Film Diffusion
% [LDG] = getLDGThinFilmMatrix(qRiemann, num_elems, num_basis_cpts, deltaX);
% spy(LDG);
% qRiemannVector = getQVector(qRiemann);
% LDGqVector = LDG*qRiemannVector;
% 
% Aq = ThinFilmDiffusionOperator(qRiemann, deltaX, qRiemann, 'extrapolation');
% AqVector = getQVector(Aq);
% norm(AqVector - LDGqVector)
% 
% Dq = ThinFilmDiffusionDiagonalOperator(ones(num_elems, num_basis_cpts), deltaX, qRiemann, 'extrapolation');
% DqVector = getQVector(Dq);
% norm(diag(LDG) - DqVector)
