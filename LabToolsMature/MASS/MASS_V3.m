% This code is created by Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan, Chetan Kumar Gupta and Eamonn Keogh.
% The overall time complexity of the code is O(n log n). The code is free to use for research purposes.
% The code may produce imaginary numbers due to numerical errors for long time series where batch processing on short segments can solve the problem.

%x is the long time series
%y is the query
%k is the size of pieces, preferably a power of two

function [dist] = MASS_V3(x, y, k)
%x is the data, y is the query
m = length(y);
n = length(x);
dist = [];

%compute y stats -- O(n)
meany = mean(y);
sigmay = std(y,1);

%compute x stats -- O(n)
meanx = movmean(x,[m-1 0]);
sigmax = movstd(x,[m-1 0],1);

%k = 4096; %assume k > m
%k = pow2(nextpow2(sqrt(n)));

y = y(end:-1:1);%Reverse the query
y(m+1:k) = 0; %aappend zeros

for j = 1:k-m+1:n-k+1
   
   
%The main trick of getting dot products in O(n log n) time
X = fft(x(j:j+k-1));
Y = fft(y);
Z = X.*Y;
z = ifft(Z);

d = 2*(m-(z(m:k)-m*meanx(m+j-1:j+k-1)*meany)./(sigmax(m+j-1:j+k-1)*sigmay));
dist = [dist ; sqrt(d)];

end

j = j+k-m;
k = n-j;
if k >= m
   
%The main trick of getting dot products in O(n log n) time
X = fft(x(j+1:n));

y(k+1:end)= [];

Y = fft(y);
Z = X.*Y;
z = ifft(Z);

d = 2*(m-(z(m:k)-m*meanx(j+m:n)*meany)./(sigmax(j+m:n)*sigmay));
dist = [dist ; sqrt(d)];

end