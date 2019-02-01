f = @(x) sin(2*pi*x);
exactSol = @(x) sin(2*pi*x)/(4*pi^2);
exactSolD = @(x) -cos(2*pi*x)/(2*pi);
a = 0.0;
b = 1.0;

error = [];
errorD = [];
hArray = [];
for N = [20, 40, 80, 160, 320]
    h = (b-a)/N;
    hArray = [hArray, h];
    x = @(i) a + i*h-h/2;
    leftbc = 0.0;
    rightbc = 0.0;
    tau = 1.0;
    
    e = ones(N,1);
    cl = tau/2.0 + 1.0/h;
    cc = -1.0*tau - 2.0/h;
    cr = tau/2.0 + 1.0/h;
    A = spdiags([cl*e,cc*e, cr*e],-1:1,N-1,N-1);
    rhs = (-h/2)*(f(x(1:(N-1))) + f(x(2:N)))';
    rhs(1) = rhs(1) - cl*leftbc;
    rhs(end) = rhs(end) - cr*rightbc;
    
    uhat = A\rhs;
    uhat = [leftbc; uhat; rightbc];
    
    Qf = zeros(N,1);
    Quhat = (-1.0/h)*(uhat(2:end) - uhat(1:end-1));
    Uf = (h/(2*tau))*f(x(1:N))';
    Uuhat = (1.0/2.0)*(uhat(2:end) + uhat(1:end-1));
    
    qh = Qf + Quhat;
    uh = Uf + Uuhat;
    
    %plot(x(1:N), qh, x(1:N), exactSolD(x(1:N)))
    error = [error, norm(uh - exactSol(x(1:N))')];
    errorD = [errorD,norm(qh - exactSolD(x(1:N))')];
end
hRatios = hArray(1:end-1)./hArray(2:end);
errorRatios = error(1:end-1)./error(2:end);
errorDRatios = errorD(1:end-1)./errorD(2:end);
order = log(errorRatios)./log(hRatios);
orderD = log(errorDRatios)./log(hRatios);
table(hArray', error', [0,order]', errorD', [0,orderD]')
plot(x(1:N), uh, x(1:N), exactSol(x(1:N)))
