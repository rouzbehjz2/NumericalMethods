
% Polynomial interpolation on eq. nodes
clear all; close all;
xbar = linspace(-1,1,200)'; % needed for plotting the fun
ybar = 1./(1+6*xbar.^2); % needed for plotting the fun
s = 1;
for n=3:3:39
    close all
    x =linspace(-1,1,n+1);    %equispaced interpolation pts
    y =1./(1+6*x.^2);         %Runge fct at the interp. pts
    String = ['The conditioning of the Vandermonde matrix with n=' ,num2str(n+1),' equispaced data is:']
    V = vander(x);
    % V c = y 
    % c = V\y;
    Cond_number(s) = cond(V,1); 
    
    Cond_number(s)

% Compute the coefficients via polyfit     
    c = polyfit(x,y,n);
    p = polyval(c,xbar);

% Alternative 
% Compute the coefficients via the Vandermonde matrix    
%     c = V\y';
%     p = polyval(c,xbar);

% Compute the maximum error
    max_err(s) = max(abs(p-ybar));
    s = s+1;
    
    figure
    hold on
    box on
    axis square    
    plot(xbar,ybar,'g--', x,y,'bo',xbar,p,'r-');
    grid on
    title(['n=' num2str(n+1) ' data'])
    legend('Original function','Interpolation points',...
        'Interpolation on equispaced points'), pause(2.5)
end
%%
figure 
semilogy(3:3:39,Cond_number,'-*m')
title(['Condition number of the Vandermonde matrix (Equispaced data)'])
hold on 
box on 
axis square 
grid on 

%%
figure 
semilogy(3:3:39,max_err,'-*m')
title(['Max interpolation error (Equispaced data)'])
title(['Maximum Error'])
xlabel('n')
ylabel('Maximum error')
hold on 
box on 
axis square 
grid on 
