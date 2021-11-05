% Input:
%   data              is the raw time series. 
%   N                 is the length of sliding window (use the length of the raw time series
%                     instead if you don't want to have sliding windows)
%   n                 is the number of symbols in the low dimensional approximation of the sub sequence.

function [PAA] = newpaa(data, numCoeff)

N = length(data);
n = numCoeff;
win_size = floor(N/n);                         

for i = 1 : length(data) - (N -1)                                       

    sub_section = data(i:i + N -1); 

    % N is not dividable by n
    if (N/n - floor(N/n))                               
        temp = zeros(n, N);
        for j = 1 : n
            temp(j, :) = sub_section;
        end
        expanded_sub_section = reshape(temp, 1, N*n);
        PAA = [mean(reshape(expanded_sub_section, N, n))];
    % N is dividable by n
    else                                                
        PAA = [mean(reshape(sub_section,win_size,n))];
    end
    PAA = PAA';
end