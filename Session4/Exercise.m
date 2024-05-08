% Exercise 4: Eigenvalues and Singular Values

clc
disp('')
disp('exercise_4_1')
disp('')

clear all
close all
m_max = 100;
tol = 1.0e-10;
z = [1;2;3];
A_1 = [1 2 0; 1 0 0; 0 1 0];
[lambda_1,w1,m1] = powerIt(A_1,z,tol,m_max)
plot(1:m1+1,lambda_1,'b+-','linewidth',2)
lambda_max_1=lambda_1(end)
eigenvalues_A_1_eig = eig(A_1)
% The method converges quite quickly because the speed
% depends on the ratio |lambda_2|/|lambda_1|=1/2 and, therefore,
% reaches the tolerance tol in 35 iterations
pause
%%%%%
A_2 = [0.1 3.8 0; 1 0 0; 0 1 0];
[lambda_2,w2,m2] = powerIt(A_2,z,tol,m_max)
plot(1:m2+1,lambda_2,'b+-','linewidth',2)
lambda_max_2 = lambda_2(end)
eigenvalues_A_2_eig = eig(A_2)
% The method converges slowly because the speed
% depends on the ratio |lambda_2|/|lambda_1|=1.9/2=0.95
% and, therefore, 100 iterations are not sufficient to
% reach the specified tolerance
pause
%%%%%
A_3 = [0 -1 0; 1 0 0; 0 1 0];
[lambda_3,w3,m3] = powerIt(A_3,z,tol,m_max)
plot(1:m3+1,lambda_3,'b+-','linewidth',2)
lambda_max_3 = lambda_3(end)
eigenvalues_A_3_eig = eig(A_3)
% The method does not converge because A has two complex
% conjugate eigenvalues i and -i with maximum modulus

disp('END OF EXERCISE')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
disp('')
disp('exercise_4_2')
disp('')

clear all
close all
m_max = 100;
tol = 1.0e-10;
z = ones(3,1);
p = 0.5;
A_1 = [1 -2 0; 0 2 0; 1 1 3];
[lambda_1,w1,m1] = InvPowerIt(A_1,p,z,tol,m_max)
plot(1:m1+1,lambda_1,'b+-','linewidth',2)
eigenvalue_p = lambda_1(end)
all_eigenvalues_A_1_eigs = eigs(A_1,3)
eigenvalues_A_1_eigs_p = eigs(A_1,1,p)
% The inverse power method converges to the eigenvalue 1,
% which is the eigenvalue of A closest to p.
pause

A_2 = [0.5 -2 0; 0 2 0; 1 1 3];
[lambda_2,w2,m2] = InvPowerIt(A_2,p,z,tol,m_max)
plot(1:m2+1,lambda_2,'b+-','linewidth',2)
eigenvalue_p = lambda_2(end)
all_eigenvalues_A_2_eigs = eigs(A_2,3)
eigenvalues_A_2_eigs_p = eigs(A_2,1,0.49)
% The inverse power method does not converge because p=0.5
% is an eigenvalue of A, and thus (A-pI) is not invertible!!!
% It would be necessary to check for the possibility that p might
% already be an eigenvalue of A; in that case, it doesn't make
% sense to apply the method.
pause

A_3 = [0 -2 0; 0 1 0; 1 1 3];
[lambda_3,w3,m3] = InvPowerIt(A_3,p,z,tol,m_max)
plot(1:m2+1,lambda_3,'b+-','linewidth',2)
eigenvalue_p = lambda_3(end)
all_eigenvalues_A_3_eigs = eigs(A_3,3)
eigenvalues_A_3_eigs_p = eigs(A_3,2,p)
% The inverse power method does not converge because p=0.5
% is equidistant from the eigenvalues 0 and 1 of matrix A_3,
% and thus there exist two distinct eigenvalues (real and
% opposite in sign) with maximum modulus for the matrix (A-pI)^(-1)!

disp('END OF EXERCISE')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
disp('')
disp('exercise_4_3')
disp('')

clear all
close all
n = 10;
A = hilb(n);
m_max = 100;
tol = 1.0e-14;
[d,m] = qr_base(A,tol,m_max);
err = abs(sort(eig(A))-sort(d))

% The maximum absolute error associated with the eigenvalues
% determined using the QR method is 1.3323e-15.

B1 = [0 0 2; 1 0 1; 0 1 1];
n = size(B1,1);
b = eig(B1);
% The matrix B has one real eigenvalue and two complex conjugate eigenvalues.
% Convergence of the QR method in its basic form is not guaranteed. Nonetheless,
% we perform 100 iterations of the method.
A = B1;
for m = 1:m_max
[Q,R] = qr(A);
A = R*Q;
end
A
% The structure of the matrix A obtained at the 100th iteration
% is nearly triangular. We verify that the eigenvalues deduced
% from the final matrix A are indeed approximations of the
% eigenvalues of B.
a = zeros(n,1);
a(1) = A(1,1);
a(2:3) = eig(A(2:3,2:3));
err = abs(sort(a)-sort(b))
% Thus, the basic form of the QR method has provided an approximation
% of the eigenvalues of matrix B.

disp('END OF EXERCISE')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
disp('')
disp('exercise_4_4')
disp('')

clear all
close all
format long e
for n = 5:5:100
A = -ones(n);
A = triu(A,1)+diag(ones(n,1));
[U,S,V] = svd(A);
determinant = det(A)
rank_A = rank(A)
s = diag(S);
pause
end
% The determinant of matrix A is 1 (product of the diagonal elements)
% and therefore, the effective rank of the matrix is n.
% The numerical rank (which is calculated up to a tolerance) is, however, n-1,
% starting from a certain order n of the matrix onwards!
% This result should be interpreted as follows:
% as n increases, the matrix approaches a singular matrix more and more.
% In fact, the last singular value s(n) decreases as n increases,
% confirming that the given matrix becomes closer to the set of matrices
% with rank n-1 as n increases.
% It should be noted that this behavior cannot be deduced from the value of
% the determinant, which remains constantly equal to 1, but from the numerical rank,
% i.e., the singular values of the matrix.

disp('END OF EXERCISE')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
disp('')
disp('exercise_4_5')
disp('')

clear all
close all
format short e
A = [3 -2 1 2; -1 0 2 1; 0 5 -6 -1; 1 1 -1 1; 1 -1 -1 -1; 8 -1 -5 2];
b = [1; -3; 7; 0; -6; 2];
r = rank(A); % The matrix has rank 3
rank_Ab = rank([A b]); % The augmented matrix has rank 4
% The system does not have a classical solution, and
% we compute the least squares solution.
% Since the matrix A does not have full rank 4,
% there exist infinitely many least squares solutions,
% but only one of them has the minimum 2-norm. This solution
% can be obtained using the singular value decomposition of A.
[U,S,V] = svd(A);
singular_values = diag(S);
ystar = zeros(4,1);
ystar(1:r) = (U(:,1:r)'*b)./singular_values(1:r);
xstar = V*ystar;
% The computed solution can also be obtained
% simply with the following instruction.
x = pinv(A)*b;
err = norm(xstar-x)
%
% Note that all the least squares solutions of the system
% are given by xstar+z, where z belongs to the kernel of A.
% The kernel of A is generated by the fourth column vector of V,
% and therefore, the solutions are of the form xstar+cV(:,4),
% where c is a constant.

disp('END OF EXERCISE')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

disp('*******************************************')
disp('exercise_4_6_Quiz')
disp('**************************************')

disp('The condition number in 2-norm is obtained from the ratio between the maximum and minimum singular values, which is 20')

disp('END OF EXERCISE')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

disp('*******************************************')
disp('exercise_4_7_Quiz')
disp('**************************************')

A=[1 5 0 3 9; 7 8 4 0 1; 2 5 3 9 0; 1 -1 2 1 1; 7 3 -2 0 1];

[Q,R]=qr(A);
for i=1:100
A=R*Q;
[Q,R]=qr(A);
end
A
% It can be observed that the matrix has a diagonal 2x2 block at positions 3:4

disp('END OF EXERCISE')
pause


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDITIONAL EXERCISES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
disp('')
disp('****exercise_4_1_ADD')
disp('')

clear all
close all
format long e
n = 100;

% To generate A1, we proceed as follows:
% 1) Generate the vector p whose elements are n,n-1,...,2,1
% 2) Repeat the vector p n times using the command p(ones(n,1), :)
% to generate a matrix whose n rows are all equal to the vector p
% 3) Generate a matrix with all elements equal to zero except the
% ones on the subdiagonal, which are equal to 1, using the command
% diag(ones(n-1,1), -1)
% 4) Extract the elements of A starting from the subdiagonal and set
% the ones below the subdiagonal to 0 using triu(A,-1)

p = n:-1:1;
A1 = triu( p( ones(n,1), :) - diag( ones(n-1,1), -1), -1 );

% Compute the eigenvalues of A1
lambda1 = eig(A1);
% Plot the eigenvalues
plot(real(lambda1),imag(lambda1),'r*','markersize',6)
hold on
% Generate the perturbed matrix A1p by perturbing the elements of the last row of A1
A1p = A1;
A1p(n,:) = A1p(n,:)+1.0e-10;
% Compute the eigenvalues of A1p
lambda1p = eig(A1p);
% Plot the eigenvalues
plot(real(lambda1p),imag(lambda1p),'ko','markersize',6)
pause

% From the plot, it can be deduced that the problem is ill-conditioned
% since a small perturbation in the elements of the matrix does not correspond
% to a perturbation in the eigenvalues of the same order of magnitude

% The ill-conditioning is indicated by the condeig command:
% the components of the output vector, which represent the condition
% number of each eigenvalue, are much greater than 1
condition_number_A1_eigenvalues = condeig(A1)
hold off
close
pause

% Repeat the same instructions as before for matrix A2
A2 = triu(A1)+triu(A1,1)';
lambda2 = eig(A2);
plot(real(lambda2),imag(lambda2),'r*','markersize',6)
hold on
A2p = A2;
A2p(n,:) = A2p(n,:)+1.0e-10;
lambda2p = eig(A2p);
plot(real(lambda2p),imag(lambda2p),'ko','markersize',6)
pause
condition_number_A2_eigenvalues = condeig(A2)
err_matrix_A2 = norm(A2-A2p)
err_eig_A2 = abs(sort(lambda2)-sort(lambda2p))
% In this case, the problem is well-conditioned (the matrix is symmetric)
% and the eigenvalues of the perturbed matrix are affected by an error
% of the same order of magnitude (even smaller) as the error associated
% with the perturbed matrix

disp('END OF EXERCISE')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
disp('')
disp('****exercise_4_2_Add')
disp('')

clear all
close all
z = (1:3)';
A = [0.1 3.8 0; 1 0 0; 0 1 0];
for m_max = 500:300:1100
[lambda,w,m] =powerit_no_norm(A,z,m_max);
w
pause
end
% If the iterate vectors are not normalized, there are overflow problems!

disp('END OF EXERCISE')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('')
disp('****exercise_4_3_Add')
disp('')

clear all
close all
clc
A = imread('puppy.jpg');
% The imread command reads the image saved in the file puppy.jpg
% and returns the matrix A, where each element A(i,j) contains
% the color code of the pixel (i,j) (pixel is derived from the contraction of the phrase
% "picture element" and represents the smallest element that constitutes an image).
% A is an M-by-N matrix if the image is grayscale;
% A is an M-by-N "three-plane" matrix if the image is in color
% (the first plane is for red, the second for green, and the third for blue).
%
B = rgb2gray(A);
% The rgb2gray command converts the color image to grayscale.
% RGB is an "additive" color model, which means it is a system
% based on three primary colors
% (which should not be confused with primary colors)
% that are red, green, and blue, and three subtractive colors,
% which are yellow, magenta, and cyan.
% RGB stands for Red, Green, and Blue, the names of the additive colors in English.
%
size(B)
figure
imshow(B) % to display the image
pause
B = double(B); % convert B to double precision
[U,S,V] = svd(B);
for n = 10:20:70
An = U(:,1:n)*S(1:n,1:n)*V(:,1:n)'; % rank-n approximation of matrix B
figure
Bn = uint8(An); % convert to integer
imshow(Bn)
title(['n=',int2str(n)])
end
disp('END OF EXERCISE')
