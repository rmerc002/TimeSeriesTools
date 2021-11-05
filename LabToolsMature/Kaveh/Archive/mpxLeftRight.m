function [rMatrixProfile, lMatrixProfile, rMatrixProfileIdx, lMatrixProfileIdx] = mpxLeftRight(timeSeries, minlag, subseqLen)



% Code and skew symmetric update formulas by Kaveh Kamgar. They originated
% as a modification to earlier work by Yan Zhu.
% GUI and top k motif critera are set up to match the prior work by Michael Yeh as close as possible.
%
% Additional References
% Yan Zhu, et al, Matrix Profile II: Exploiting a Novel Algorithm and GPUs to break the one Hundred Million Barrier for Time Series Motifs and Join
% Zachary Zimmerman, et al, Scaling Time Series Motif Discovery with GPUs: Breaking the Quintillion Pairwise Comparisons a Day Barrier. (pending review)
% Philippe Pebay, et al, Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments
% Takeshi Ogita, et al, Accurate Sum and Dot Product
%

n = length(timeSeries);

% difference equations have 0 as their first entry here to simplify index
% calculations slightly. Alternatively, it's also possible to swap this to the last element
% and reorder the comparison step (or omit on the last step). This is a
% special case when comparing a single time series to itself. The more general
% case with time series A,B can be computed using difference equations for
% each time series.

if nargin ~= 3
    error('incorrect number of input arguments');
elseif ~isvector(timeSeries)
    error('first argument must be a 1D vector');
elseif ~(isfinite(subseqLen) && floor(subseqLen) == subseqLen) || (subseqLen < 2) || (subseqLen > length(timeSeries))
    error('subsequence length must be an integer value between 2 and the length of the timeSeries');
end

transposed_ = isrow(timeSeries);
if transposed_
    timeSeries = transpose(timeSeries);
end

nanIdx = isnan(timeSeries);
nanSubseqs = movsum(~isfinite(timeSeries), [0 subseqLen-1], 'Endpoints', 'discard') ~= 0;
rMatrixProfile = repmat(-1, n - subseqLen + 1, 1);
rMatrixProfile(nanSubseqs) = NaN;
lMatrixProfile = repmat(-1, n - subseqLen + 1, 1);
lMatrixProfile(nanSubseqs) = NaN;
timeSeries(nanIdx) = 0;
mu = moving_mean(timeSeries, subseqLen);
invnorm = fastinvn(timeSeries, mu, subseqLen);
invnorm(nanSubseqs) = NaN;


% The update equations diff_f and diff_g are slightly difficult to derive directly, but they arise naturally taking the difference
% between cov(x(i + 1 : i + subseqLen), x(j + 1 : j + subseqLen)) and
% cov(x(i : i + subseqLen - 1), x(j : j + subseqLen - 1))
% from Pebay et al, then applying the identity
% ac - bd = (1/2) * ((a + b) * (b - d) + (a - b) * (b + d))
% to the resulting components.

% Covariance updates which involve a high level of numerical cancellation
% then correspond to very very small changes in covariance and seem to
% rarely ever become a problem in practice, particularly as motif discovery isn't as
% sensitive to perturbations as the algorithms used here to compute cross correlation and euclidean distance.
df = [0; (1/2)*(timeSeries(1 + subseqLen : n) - timeSeries(1 : n - subseqLen))];
dg = [0; (timeSeries(1 + subseqLen : n) - mu(2 : n - subseqLen + 1)) + (timeSeries(1 : n - subseqLen) - mu(1 : n - subseqLen))];

lMatrixProfileIdx = NaN(n - subseqLen + 1, 1);
rMatrixProfileIdx = NaN(n - subseqLen + 1, 1);
comparesTotal = (n - subseqLen + 1) * (n - subseqLen - minlag + 1) / 2;  % <-- each update actually updates 2
updateRate = comparesTotal/100;
comparesSoFar = 0;
counter = 0;

% The terms row and diagonal here refer to a hankel matrix representation of a time series
% This uses scaled cross correlation (scaled to allow the use of movstd) as an intermediate quantity for performance reasons.
% It is later reduced to z-normalized euclidean distance.
for diag = minlag + 1 : n - subseqLen + 1
    if comparesTotal > 32768 && counter >= updateRate
        comparesSoFar = comparesSoFar + counter;
%         fprintf('we are approximately %d percent complete\n', floor(100 * comparesSoFar/comparesTotal));
        counter = 0;
    end
    cov_ = (sum((timeSeries(diag : diag + subseqLen - 1) - mu(diag)) .* (timeSeries(1 : subseqLen) - mu(1))));
    for row = 1 : n - subseqLen - diag + 2
        cov_ = cov_ + df(row) * dg(row + diag - 1) + df(row + diag - 1) * dg(row);
        corr_ = cov_ * invnorm(row) * invnorm(row + diag - 1);
        if corr_ > rMatrixProfile(row)
            rMatrixProfile(row) = corr_;
            rMatrixProfileIdx(row) = row + diag - 1;
        end
        if corr_ > lMatrixProfile(row + diag - 1)
            lMatrixProfile(row + diag - 1) = corr_;
            lMatrixProfileIdx(row + diag - 1) = row;
        end
    end
    counter = counter + (n - subseqLen - diag + 2);
    
end

% This is the correct ordering. findMotifsDiscords uses a correlation based
% updates to avoid it being problematically slow if called in a loop.
%[motifsIdx, mpAugmented] = findMotifs(timeSeries, mu, invnorm, matrixProfile, matrixProfileIdx, subseqLen, 3, 10, minlag, 2);
%[discordsIdx] = findDiscords(mpAugmented, 3, minlag);
%timeSeries(nanIdx) = NaN;
rMatrixProfile = sqrt(2 * subseqLen * (1 - min(1, rMatrixProfile, 'includenan')));
lMatrixProfile = sqrt(2 * subseqLen * (1 - min(1, lMatrixProfile, 'includenan')));


if transposed_  % matches the profile and profile index but not the motif or discord index to the input format
    rMatrixProfile = transpose(rMatrixProfile);
    rMatrixProfileIdx = transpose(rMatrixProfileIdx);
    lMatrixProfile = transpose(lMatrixProfile);
    lMatrixProfileIdx = transpose(lMatrixProfileIdx);
end

end

function [invn] = fastinvn(ts, mu, sublen)
% This is a simple variation on Welford's method.
% This version still results in some minor cancellation and could be
% improved. It isn't prone to anything disastrous.
invn = zeros(length(mu), 1);
invn(1) = sum((ts(1 : sublen) - mu(1)).^2);
for i = 2:length(invn)
    invn(i) = invn(i - 1) + ((ts(i - 1) - mu(i - 1)) + (ts(i + sublen - 1) - mu(i))) * (ts(i + sublen - 1) - ts(i - 1));
end
invn = 1./sqrt(invn);
end


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
