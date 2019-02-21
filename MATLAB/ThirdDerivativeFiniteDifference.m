% this shows that this finite difference is second order
f = @(x) exp(x);
f_exact_third_derivative = @(x) exp(x);
N = 1000;
error = zeros(N,1);
hArray = logspace(-5,0,N);
for i = 1:N
    h = hArray(i);
    a = 1;
    f_im1 = f(a - 3*h/2);
    f_i = f(a - h/2);
    f_ip1 = f(a + h/2);
    f_ip2 = f(a + 3*h/2);
    finite_difference = (1/h^3)*(-f_im1 + 3*f_i - 3*f_ip1 + f_ip2);
    error(i) = abs(finite_difference - f_exact_third_derivative(a));
end

loglog(hArray, error);
% shows that slope is 2
polyfit(log(hArray(500:end)),log(error(500:end)'),1)