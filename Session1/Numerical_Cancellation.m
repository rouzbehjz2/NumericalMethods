%Exercise: How many strictly positive machine numbers do exist in a machine working with base N=2, t=3 and -2 <= q <= 3 and rounding to even technique.

N=2; % only digits 0 and 1
t=3; % 3 digits of mantissa
L=-2; % lower limit for exponent
U=3; % upper limit for exponent
% 0.1a2a3*N^(q) is the form of the number (different from 0, as the exercise is asking for strictly positive numbers)
% Mantissa can be seen as the sum of digits times the exponent indicating the position: for example
%number= [1,1,1]*(2.^[-1,-2,-3])';

%the only possibilities for mantissa are [110],[111],[100][101]
for i=0:1
    for q=-2:3
        matrix_number(i+1,q+3)=(2^q)*([1,0,i]*(2.^[-1,-2,-3])');
        matrix_number(i+3,q+3)=(2^q)*([1,1,i]*(2.^[-1,-2,-3])');
    end
end
matrix_number
%As can be seen by the output, numbers are 24. You can also use the command isequal to check that are all different. 

