function [q] = reimannQ(nBasisCpts, nCells, qLeft, qRight)
    q = zeros(nBasisCpts, nCells);
    q(1,1:floor(nCells/2)) = qLeft;
    q(1,floor(nCells/2))+1:end) = qRight;
end
