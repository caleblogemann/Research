function [A] = getAMatrix(nBasisCpts)
    SMat = [0, 0, 0, 0, 0;
    2*sqrt(3), 0, 0, 0, 0;
    0, 2*sqrt(3)*sqrt(5), 0, 0, 0;
    2*sqrt(7), 0, 2*sqrt(5)*sqrt(7), 0, 0;
    0, 6*sqrt(3), 0, 6*sqrt(7), 0];
    A = sparse(SMat[1:nBasisCpts,1:nBasisCpts]);
end

%{{{{0.0, 0.0, 0.0, 0.0, 0.0}},
%{{2.0*CONSTANTS::SQ3, 0.0, 0.0, 0.0, 0.0}},
%{{0.0, 2.0*CONSTANTS::SQ3*CONSTANTS::SQ5, 0.0, 0.0, 0.0}},
%{{2.0*CONSTANTS::SQ7, 0.0, 2.0*CONSTANTS::SQ5*CONSTANTS::SQ7, 0.0, 0.0}},
%{{0.0, 6.0*constants::sq3, 0.0, 6.0*constants::sq7, 0.0}}}}; 
