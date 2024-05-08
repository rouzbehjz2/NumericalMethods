function [d,m] = qr_base(A,toll,m_max)
%convergence is guaranteed only if A has eigenvalues with distinct absolute value
for m = 1:m_max
    [Q,R] = qr(A);
    A = R*Q;
    if norm(tril(A,-2),inf) <= toll
       break
    end
end
d = diag(A);
