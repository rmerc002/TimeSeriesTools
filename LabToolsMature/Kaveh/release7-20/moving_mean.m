
function [ res ] = moving_mean(a,w)
% moving mean over sequence a with window length w
% based on Ogita et. al, Accurate Sum and Dot Product

% This moving mean function is used to reduce accumulated rounding error (in spite of the fact that motif
% discovery isn't that sensitive to them).
% Using a typical summation algorithm to compute the mean for each window independently, the
% resulting errors attributable to accumulated roundoff form a not
% necessarily stationary time series. By our observations, this is typically the largest source
% of accumulated rounding error in the use of iterative accumulation algorithms to compute cross covariance.

% Using compensated arithmetic to compute the mean, it takes less time than
% computing each window separately, and the results of cross covariance re
% much closer to a naives evaluation of dot(x - mu_x, y - mu_y), which is much less sensitive to minor perturbations in the mean than
% any online method for computing covariance. The actual error bounds are still linear with respect to inputs and conditioning,
% but the worst case suggested by forward error analysis is much nicer when we can assume that the computed mean is very close
% to the closest floating point approximation of the mean.


res = zeros(length(a) - w + 1, 1);
p = a(1);
s = 0;

for i = 2 : w
    x = p + a(i);
    z = x - p;
    s = s + ((p - (x - z)) + (a(i) - z));
    p = x;
end

res(1) = p + s;

for i = w + 1 : length(a)
    x = p - a(i - w);
    z = x - p;
    s = s + ((p - (x - z)) - (a(i - w) + z));
    p = x;
    
    x = p + a(i);
    z = x - p;
    s = s + ((p - (x - z)) + (a(i) - z));
    p = x;
    
    res(i - w + 1) = p + s;
end

res = res ./ w;

end
