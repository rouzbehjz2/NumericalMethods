function fx = val_f(x)
%function used in exercise 0_2 ADD
%Input
%x: vector of points at which evaluate f
%Output
%fx: vector which contains values of f at points of vector x
%
n = length(x);
fx = zeros(1,n);

% evaluate f according to x
for i = 1:n
    if x(i) < 0
        fx(i) = -2*x(i);
    elseif x(i) == 0
        fx(i) = 0;
    elseif x(i) > 0
        fx(i) = 2*x(i);
    end
end
