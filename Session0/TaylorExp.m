function [v,i] = taylor_exp(x,toll)
% function used in exercise 0_1 Add
%Input
%x: given point at which evaluate the function
%toll: tolerance
%Output
%v: value of the polynomial at x   
%i: degree of the a polynomial
%term: term of the Taylor series 
v = 0;
i = 0;
term = 1;
while term >= toll
    v = v+term;
    i = i+1;
    term=x^i/factorial(i);
end
