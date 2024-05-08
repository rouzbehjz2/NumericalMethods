%Session 2: Data and Function approximation

clear all
close all
clc
disp('***********************************************')
disp('*****************EXERCISE_2_1*****************')
disp('***********************************************')

a = 0;
b = pi;
f = @(x) sin(x); 
%a = -2*pi;
%b = 2*pi;
%f = @(x) 1./(1+x.^2); 
z = linspace(a,b);
fz = f(z);
for n = 5:5:15
   x = linspace(a,b,n+1); %equispaced nodes
   y = f(x);
   c1 = polyfit(x,y,n);
   p1 = polyval(c1,z);
   err1 = norm(p1-fz,inf)
   err1_100 = abs(p1-fz);
   t = -cos(((1:n+1)-1)*pi/n); 
   x = (b-a)/2*t+(b+a)/2; %Chebyshev-Lobatto nodes scaled and shifted
   y = f(x);
   c2 = polyfit(x,y,n);
   p2 = polyval(c2,z);
   err2 = norm(p2-fz,inf)
   err2_100 = abs(p2-fz);
   figure(1)
   plot(z,p1,'b',z,p2,'r',z,fz,'g',x,y,'ko','linewidth',3)
   legend('equispaced','Chebyshev-Lobatto','sin(x)','data')
   figure(2)
   plot(z,err1_100,'b',z,err2_100,'r','linewidth',3)
   legend('equispaced','Chebyshev-Lobatto')
   pause
end
disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear all
clc
disp('***********************************************')
disp('*****************EXERCISE_2_2*****************')
disp('***********************************************')
f = @(x) (1-x.^2).^(5/2);
fd = @(x) (5/2)*(1-x.^2).^(5/2-1).*(-2*x);
f0 = fd(-1);
fn = fd(1);
z = linspace(-1,1);
fz = f(z);
for k = 2:5
    n = 2^k;
    x = -1+2*(0:n)/n;
    y = f(x);
    s = spline(x,y,z);
    s1 = spline(x,[f0 y fn],z);
    figure(1)
    plot(x,y,'ko',z,fz,'r',z,s,'b',z,s1,'g','linewidth',3)
    legend('data','f(x)','not-a-knot spline ','constrained spline')
    pause
    figure(2)
    semilogy(z,abs(s-fz),'b',z,abs(s1-fz),'g','linewidth',3)
    legend('error not-a-knot spline ','error constrained spline ')
    err = norm(fz-s,inf)
    err1 = norm(fz-s1,inf)
    pause
end    
% constrained spline gives a more accurate approximation 
% with respect to a not-a-knot spline because, unlike
% the latter, in addition to the interpolation conditions, satisfies
% two further conditions linking it to the function f

disp('**********************************************FINE EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear all
clc
disp('***********************************************')
disp('*****************EXERCISE_2_3_QUIZ************')
disp('***********************************************')
fun=@(x) cos(x)./(1+3*x.^2);
a=-1;
b=1;
n=7;
xe=linspace(a,b,n+1);
ii=1:n+1;
z=-cos((2*ii-1)*pi/(2*(n+1)));
xc=0.5*(a+b)+0.5*(b-a)*z;

pe=polyfit(xe,fun(xe),n);

pc=polyfit(xc,fun(xc),n);

xvis=linspace(a,b,1000); % I use a high number of nodes to display the functions well

% To evaluate the infinite norm of the error (i.e. I have to find the maximum in module) I graph the error in module
%e I proceed to visually identify the maximum and the value of the abscissa for which this maximum is obtained
figure(1)
plot(xvis,abs(fun(xvis)-polyval(pe,xvis)),'k')
% From the graph it can be seen that the maximum error for p and x is around +/-0.9
% from which:
Ee=abs(fun(0.9)-polyval(pe,0.9)); % Attention: observe that the infinite norm of the function fun(x)
%on the interval is 1, so the absolute error coincides with the relative error!!!

figure(2)
plot(xvis,abs(fun(xvis)-polyval(pc,xvis)),'k')

% From the graph it can be seen that the maximum error per pc occurs for x = 0
% from which:
Ec=abs(fun(0)-polyval(pc,0)); % Attention: observe that the infinite norm of the function fun(x)
%on the interval is 1, so the absolute error coincides with the relative error!!!

Ee/Ec    

disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
disp('***********************************************')
disp('*****************EXERCISE_2_4_QUIZ************')
disp('***********************************************')
fun=@(x) cos(x).*log(x);

a=1;
b=5;
n=3;
x=linspace(a,b,n+1);

p=polyfit(x,fun(x),n);

c=[1/4 1/3 1/2 1]; % vector of the coefficients of the monomials that corresponds to the operation of integration of the polynomial: i.e. the integral of x^3=1/4 x^4, and so on, obtaining [1/4 1/3 1/2 1]

% The c.*p operation calculates the integral. The concatenated zero at the bottom represents the monomial of degree zero - the constant, which is arbitrarily set to zero - and obviously has no effect on the computation of the definite integral. However, it must be inserted because in this way the coefficients are interpreted by the polyval as coefficients of monomials of degree from 1 to 4.

I=polyval([c.*p,0],5)-polyval([c.*p,0],1)
% Finally, I use polyval to evaluate the integrated polynomial in the extrema of integration
disp('**********************************************END EXERCISE**********************************************')
pause


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADDITIONAL EXERCISES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc
format short e

disp('***********************************************')
disp('*****************EXERCISE_2_1_ADD*************')
disp('***********************************************')

a = -5;
b = 5;
f = @(x) 1./(1+x.^2); 
z = linspace(a,b);
fz = f(z);
for n = 5:4:13
   x=linspace(a,b,n+1); %equispaced nodes
   %t = -cos((2*(1:n+1)-1)*pi/(2*(n+1))); x = (b-a)/2*t+(b+a)/2; % Chebyshev nodes scaled and shifted
   y = f(x);
   c = polyfit(x,y,n);
   p = polyval(c,z);
   plot(z,p,'b',z,fz,'r',x,y,'ko','linewidth',3)
   pause
end

disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
disp('***********************************************')
disp('*****************EXERCISE_2_2_ADD*************')
disp('***********************************************')

a = -5;
b = 5;
f = @(x) 1./(1+x.^2); 
z = linspace(a,b);
fz = f(z);
for n = 6:4:14
    x = linspace(a,b,n);
    y = f(x);
    s = spline(x,y,z);
    plot(x,y,'ko',z,fz,'r',z,s,'b','linewidth',3)
    err = norm(s-fz,inf)
    pause
end    

disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear all
clc
disp('***********************************************')
disp('*****************EXERCISE_2_3_ADD*************')
disp('***********************************************')

x = 1:20;
y = [243 209 181 179 180 166 163 157 187 192 138 95 56 32 21 12 11 61 146 186];
z = linspace(1,20,1000);
s1 = interp1(x,y,z);
plot(z,s1,'b',x,y,'ko','linewidth',3)
hold on
s11 = interp1(x,y,[2.5 19.5])
plot([2.5 19.5],s11,'go','linewidth',3)
grid on
disp('**********************************************END EXERCISE**********************************************')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
disp('***********************************************')
disp('*****************EXERCISE_2_4_ADD_Quiz********')
disp('***********************************************')

x = [0 1 2];
y = [1 0 0];
l1=polyfit(x,y,2);
polyval(l1,3)

disp('**********************************************END EXERCISE**********************************************')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
