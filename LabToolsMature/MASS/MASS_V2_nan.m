function [dist] = MASS_V2_nan(x, y)

nanIndices = ones(length(x),1);
nanIndices(isnan(x)) = nan;

mm = length(y);

if any(isnan(x))
    x(isnan(x)) = mean(x,'omitnan');
end
%x is the data, y is the query
m = length(y);
n = length(x);

%compute y stats -- O(n)
meany = mean(y);
sigmay = std(y,1);

%compute x stats -- O(n)
meanx = movmean(x,[m-1 0]);
sigmax = movstd(x,[m-1 0],1);

y = y(end:-1:1);%Reverse the query
y(m+1:n) = 0; %aappend zeros

%The main trick of getting dot products in O(n log n) time
X = fft(x);
Y = fft(y);
Z = X.*Y;
z = ifft(Z);

dist = 2*(m-(z(m:n)-m*meanx(m:n)*meany)./(sigmax(m:n)*sigmay));
dist = sqrt(dist);
dist = dist./sqrt(length(y));

%added 2021-09-01 by Ryan to match standard mp euclidean distances
dist = real(dist);
dist = dist*sqrt(length(x));
dist = applySubsequenceNan(nanIndices, dist, mm);
end

