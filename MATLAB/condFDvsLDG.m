exactSolution = @(x, t) 0.2*exp(-10*t).*exp(-300.0*(x - 0.5).^2) + 0.1;

a = 0.0;
b = 1.0;

diffusivity = 1.0;
bc = 'periodic';

nCellsArray = [20, 40, 80];
condFD = [];
condLDG = [];

for nCells = nCellsArray
    
    deltaX = (b-a)/nCells;
    x = (a + deltaX/2):deltaX:(b-deltaX/2);
    
    q1 = projectQ(@(x) exactSolution(x, 0), 1, nCells, a, b);
    q2 = projectQ(@(x) exactSolution(x, 0), 2, nCells, a, b);
    FDMatrix = getFDThinFilmMatrix(q1, nCells, deltaX, diffusivity, bc);
    LDGMatrix = getLDGThinFilmMatrix(q2, nCells, nBasisCpts, deltaX, diffusivity, bc);
    
    condFD = [condFD, cond(full(FDMatrix))];
    condLDG = [condLDG, cond(full(LDGMatrix))];
    
end
table(condFD', condLDG')