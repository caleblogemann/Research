a = 0.0;
b = 1.0;
nCells = 6;
qleft = .3323;
qright = 0.1;
discontinuity = 0.5;
reimannFunction = @(x) qleft*(x <= discontinuity) + qright*(x > discontinuity);
q2 = projectQ(reimannFunction, 2, nCells, a, b);
A = getFDThinFilmMatrix(q2, nCells, deltaX, 1.0, 'periodic');
q_FD = getQFD(q2);
A*q_FD
