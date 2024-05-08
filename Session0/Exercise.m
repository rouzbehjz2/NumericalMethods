% session 0: vector and matrix manipulation, graphs of
% functions and MATLAB language

clear all
clc
disp('***********************************************')
disp('*****************exercise_0_1*****************')
disp('***********************************************')

x = [1:-0.1:0]
x([1 4 3]) 
disp('extracts the components 1, 4 and 3 of the vector x')
pause
%
x = [1:-0.1:0]
x([1:2:7 10])=zeros(1,5) 
disp('sets the components 1, 3, 5, 7 and 10 of x equal to 0') 
pause
%
x = [1:-0.1:0]
x([1 2 5])=[0.5*ones(1,2) -0.3] 
disp('sets the components 1, 2 and 5 of x equal to 0.5, 0.5 and -0.3, respectively')
pause
%
x = [1:-0.1:0]
y = x(end:-1:1) 
disp('defines a vector y, whose components are those of vector x ordered from last to first')

disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
disp('***********************************************')
disp('*****************exercise_0_2*****************')
disp('***********************************************')

A = [1:4;5:8;9:12] % A is a 3x4 matrix
size(A) 
disp('gives the number of rows and the number of columns of matrix A')
pause

A(1:2,4)
disp('extracts the submatrix of A formed by the elements belonging to the first two rows and the fourth column')
pause

A(3,2) = A(1,1) 
disp('sets the element in position (3,2) of the matrix equal to the element in position (1,1)')
pause

A(1:2,4) = zeros(2,1) 
disp('sets the elements of the submatrix of A, formed by the elements belonging to the first two rows and the fourth column, equal to zero')
pause

A = [1:4;5:8;9:12]
A(2,:) = A(2,:)-A(2,1)/A(1,1)*A(1,:) 
disp('redefines the second-row elements of A as a linear combination of the second-row and first-row elements:')
disp('precisely subtracts from the elements of the second row the elements of the first row multiplied by the factor A(2,1)/A(1,1)')

disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
disp('***********************************************')
disp('*****************exercise_0_3*****************')
disp('***********************************************')

A = [1:6;5:10;9:14;15:20]
%a)
B = A(:,6:-1:1)
pause
%b)
B = A(:,2:2:end)
pause
%c)
B = A(1:2:end,:)
pause
%d)
B = A([1 4 3],[5 2])
pause
%e)
diag(A)

disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
disp('***********************************************')
disp('*****************exercise_0_4*****************')
disp('***********************************************')

D = diag(5*ones(10,1));
CS = diag(3*ones(9,1),1);
CI = diag(-1*ones(9,1),-1);
B = D+CS+CI
B([5 8],[6 9]) = 2

disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
disp('***********************************************')
disp('*****************exercise_0_5*****************')
disp('***********************************************')
x = linspace(-pi,pi);
y = sin(x);
figure
plot(x,y,'b','linewidth',2)
title('sin(x)')

pause
close
%
x = linspace(-1,1);
y = exp(x);
figure
plot(x,y,'r','linewidth',2)
title('exp(x)')

pause
close
%
x = linspace(-5,5);
y = exp(-x.^2);
figure
plot(x,y,'g','linewidth',2)
title('exp(-x^2)')

pause
close
%
x = linspace(0.001,4*pi,1000);
y = sin(x)./x;
figure
plot(x,y,'m','linewidth',2)
title('sin(x)/x')

pause
close
%
%
x = linspace(0.001,2,10000);
y = x.*sin(1./x);
figure
plot(x,y,'k','linewidth',2)
title('xsin(1/x)')


disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
disp('*****************exercise_0_6*****************')
disp('***********************************************')
x = linspace(0.1,100,10000);
y = sqrt((100*(1-0.01*x.^2).^2+0.02*x.^2)./((1-x.^2).^2+0.1*x.^2));

figure
plot(x,y,'linewidth',2)
pause

figure
loglog(x,y,'linewidth',2)

disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close
disp('*****************exercise_0_7_QUIZ************')
disp('***********************************************')
doc format
disp('Note that the format command can be used with various options to change the display format of the results. It has NO impact on the accuracy of the calculations')
disp('**********************************************EDN EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDITIONAL EXERCISES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
disp('*****************esercizio_0_1_Add*************')
disp('***********************************************')
x = 1;
fx = val_f(x)
x = linspace(-1,1,5);
fx = val_f(x)
plot(x,fx)

disp('**********************************************END EXERCISE**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
disp('*****************esercizio_0_2_Add*************')
disp('***********************************************')
x = 1;
tol = 1.0e-10;
[v,i] = taylor_exp(x,tol)
err = abs(v-exp(x))/abs(exp(x))

disp('**********************************************END EXERCISE**********************************************')
