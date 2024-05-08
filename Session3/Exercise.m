%exercise session 3: lynear systems



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
disp('***********************************************')
disp('*****************exercise_3_1*****************')
disp('***********************************************')

n = 100;
for i = 1:n  % double cycle to define A
    for j = 1:n
        A(i,j) = max([i j]);
    end
end
Kinf = cond(A,inf)
b = sum(A,2);
[L,U] = elleu(A); % computes the L and U factors without executing swaps
y = L\b;  
x_nopiv = U\y; % solves the system Ax=b without applying pivoting
err_nopiv = norm(ones(n,1)-x_nopiv,inf)/norm(ones(n,1),inf) % computes the relative error associated with the solution calculated without pivoting
[L,U,P] = lu(A); % computess L and U factors by pivoting, that means performing swapping on the rows
y = L\(P*b);
x_piv = U\y;
err_piv = norm(ones(n,1)-x_piv,inf)/norm(ones(n,1),inf) % computes the relative error associated with the solution obtained with pivoting 

% the solution obtained by the method of the elimination of
% Gauss with pivoting (i.e. obtained by factoring PA=LU)
% is more accurate than without pivoting (elleu)

disp('**********************************************END EXERCISE**********************************************')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
disp('***********************************************')
disp('*****************exercise_3_2*****************')
disp('***********************************************')

n = 100;  
for i = 1:n  % double cycle to define A
    for j = 1:n
        A(i,j) = i*max([i j]);
    end
end
determinante = det(A) % computes the determinant of the matrix to verify that it is invertible
[L,U,P] = lu(A);
inverse_c = inv(U)*inv(L)*P;  % computes the inverse of A using the factors P,L,U
inverse = inv(A);  % computes the inverse of A using MATLAB's inv command
err = norm(inverse-inverse_c,inf)/norm(inverse,inf) % computes the associated relative error in norm infinity
% the two inverse_c and inverse matrices are numerically equivalent,
% since the relative error (normally infinite) is approximately equal to or lower 
% of machine accuracy

disp('**********************************************END EXERCISE**********************************************')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

disp('***********************************************')
disp('*****************exercise_3_3*****************')
disp('***********************************************')

p = 30;
n = 100;
A = rand(n); % generates a matrix of pseudo-random numbers
b = sum(A,2);
x = zeros(n,p); % defines a matrix x of dimensions n by p, in whose i-th column the solution vector of the i-th system will be stored
tic % activates the timer to calculate the calculation time for solving linear systems using PA=LU
[L,U,P] = lu(A); % compute the PA = LU factorization once and for all!
tn = P*b;% generates the right-hand side term of the first linear system
for i = 1:p
    y = L\tn;  % solves two triangular systems at each step
    x(:,i) = U\y;  
    tn = P*x(:,i); % generates the right-hand side term of the next linear system
end
toc % stops the timer

tic % activates the timer to calculate the calculation time for solving linear systems using the built-in MATLAB command \
x(:,1) = A\b;
for i = 2:p
    x(:,i) = A\x(:,i-1);
end
toc

% for large n the calculation times of the first procedure are less than
% those of the second procedure


disp('**********************************************END EXERCISE**********************************************')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
disp('***********************************************')
disp('*****************exercise_3_4*****************')
disp('***********************************************')

n = 100;
B = diag(10*ones(n,1))+diag(5*ones(n-1,1),1)+diag(-5*ones(n-1,1),-1); 
determinante = det(B)
A = B'*B;  % since B is nonsingular, A is certainly symmetric and positive definite
Kinf = cond(A,inf)
R = chol(A); % computes the Choleski factor R such that A = R'*R
R1 = inv(R);
inversa_c = R1*R1'; % computes the inverse of A using Choleski factorization
inversa = inv(A);   % computes the inverse of A using the inv MATLAB command 
err_inv = norm(inversa-inversa_c,inf)/norm(inversa,inf)
b = sum(A,2); 
y = R'\b;  % solves the system Ax=b with Choleski factorization
x = R\y;
err_sol = norm(ones(n,1)-x,inf)/norm(ones(n,1),inf) %computes the relative error associated with the solution calculated in infinite norm
disp('**********************************************END EXERCISE**********************************************')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

disp('***********************************************')
disp('*****************exercise_3_5*****************')
disp('***********************************************')

clear all
format short e
% defines a counter i, which is incremented by one
% at each step, to store in K2(i), err1(i) and err2(i)
% the number of conditioning and errors calculated at each step

i = 1;
for n = 1000:1000:5000   % other values of n that can be used are n = 100:100:500%
    A = rand(n);
    K2(i) = cond(A);
    b = sum(A,2);
    tic  % calculates the computation time of the solution of Ax = b using PA=LU method
    [L,U,P] = lu(A);
    y = L\(P*b);
    x = U\y;
    err1(i) = norm(ones(n,1)-x)/norm(ones(n,1));  %computes the relative error associated with the solution in norm 2
    t1 = toc;
    tic  % calculates the computation time of the solution of Ax = b using A=QR method
    [Q,R] = qr(A);
    x = R\(Q'*b);
    err2(i) = norm(ones(n,1)-x)/norm(ones(n,1));   %computes the relative error associated with the solution in norm 2
    t2 = toc;
    rapporto(i) = t2/t1; % computes the ratio of the calculation time of the two procedures
    i = i+1;
end
[K2' err1' err2' rapporto']  % prints the table of computed values

% for large n the calculation times t1 of the first procedure are less than
% those t2 of the second procedure. Therefore, ratio(i) > 1 for all i
disp('**********************************************END EXERCISE**********************************************')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

disp('***********************************************')
disp('*****************exercise_3_6*****************')
disp('***********************************************')

format short e
A = [1 2 3 4;-1 0 4 1;3 5 1 0;2 -1 0 1;1 1 -1 1;2 -1 0 3];
b = [1:6]';
[m n] = size(A);
r = rank(A) % computes the rank of matrix A. N.B. the rank is maximum
[Q,R] = qr(A); % computes the QR factorization of A
Rt = R(1:n,1:n); % identifies the R tilde, non-singular matrix consisting of the first n rows and n columns
c = Q'*b;  % computes vector c
c1 = c(1:n); % identifies the vector c1 made up of the first n elements of c
cstar = Rt\c1 % solve the linear system to determine the solution in least squares sense
x = A\b % calculate the same solution with the \ built-in MATLAB command for comparison
err = norm(cstar-x)/norm(x) % computes the relative error associated with the cstar solution

disp('**********************************************END EXERCISE**********************************************')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

disp('***********************************************')
disp('*****************exercise_3_7_Quiz************')
disp('***********************************************')

A=[3 1 2 4; -1 0 7 9; 0 1 2 4; 2 4 1 1];
% Using the elleu function which implements pivoting without exchanges:
[Lns,Uns]=elleu(A);
%Or do the computations by hand!

% Use the MATLAB built-in function lu which instead uses partial pivoting
[Lpiv,Upiv,P]=lu(A);
% Multipliers are in the matrices L

Lns(3,2)*Lpiv(3,2)


disp('**********************************************END EXERCISE**********************************************')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

disp('***********************************************')
disp('*****************exercise_3_8_Quiz************')
disp('***********************************************')

disp('The matrix has no maximum rank, therefore the system admits infinite solutions, but a unique solution of minimum Euclidean norm.')


disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Additional Exercises %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
format long e
disp('***********************************************')
disp('*****************exercise_3_1_ADD*************')
disp('***********************************************')

% defines a counter i, which is incremented by one
% at each step, to store the number of in Kinf(i) and err(i).
% conditioning and the error calculated at each step
i = 1;  
ordine = 5:5:15;  % order of the system = 5,10,15 
for n = ordine
    A = hilb(n);  %  Hilbert coefficient matrix
    Kinf(i) = cond(A,inf); %computes the condition number in infinite norm
    b = sum(A,2); % defines the right-hand side such that the solution coincides with the vector with all components equal to 1
    x = A\b       % solves the system
    err(i) = norm(ones(n,1)-x,inf)/norm(ones(n,1),inf); % % computes the relative error associated to the approximation provided by the MATLAB built-in command \
    i = i+1; %increments the counter i by 1
end
[ordine' err' Kinf']  % prints the table of computed values
disp('**********************************************END EXERCISE*********************************************')
pause


clear all
close all
clc

disp('***********************************************')
disp('*****************exercise_3_2_ADD*************')
disp('***********************************************')

n = 5;
A = zeros(n);
Q = zeros(n);
A(:,1) = [4;2;1;5;-1];
A(:,2) = [1;5;2;4;0];
A(:,3) = [3;10;6;2;1];
A(:,4) = [3;1;6;2;-1];
A(:,5) = [2;-1;2;0;1];
% implements the Gram-Schmidt orthogonalization procedure
for j = 1:n
    Q(:,j) = A(:,j);
    for i = 1:j-1
        Q(:,j) = Q(:,j)-A(:,j)'*Q(:,i)*Q(:,i);
    end
    Q(:,j) = Q(:,j)/norm(Q(:,j));
end
% tests the orthogonality of Q (Q'*Q=Q*Q'=I)
errs = norm(Q'*Q-eye(5),inf)
errd = norm(Q*Q'-eye(5),inf)
disp('**********************************************END EXERCISE**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













