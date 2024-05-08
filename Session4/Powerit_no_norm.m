function [lambda,w,m] = powerIt_no_norm(A,z,m_max)
w = z;
for m = 1:m_max
    z = A*w;
    lambda = w'*z/(w'*w);
    w = z;
end
