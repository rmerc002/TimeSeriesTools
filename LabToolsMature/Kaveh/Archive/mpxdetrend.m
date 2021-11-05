function [matrixProfile, matrixProfileIdx, motifsIdx, discordsIdx] = mpxdetrend(timeSeries, minlag, subseqLen)

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

nanIdx = find(isnan(timeSeries));
nanSubseqs = find(~isfinite(movsum(timeSeries, [0 subseqLen-1], 'Endpoints', 'discard')));
matrixProfile = repmat(-subseqLen, n - subseqLen + 1, 1);
matrixProfile(nanSubseqs) = NaN;
timeSeries(nanIdx) = 0;
mu = moving_mean(timeSeries, subseqLen);
center = ceil(subseqLen / 2);
if mod(subseqLen,2) == 0
   slope = (-center + 0.5 : center - 0.5)';
else
   slope = (-center + 1 : center - 1)';
end
slope = slope ./ norm(slope);
trend = zeros(length(mu), 1);
invnorm = zeros(length(mu), 1);

for i = 1 : n - subseqLen + 1
    z = (timeSeries(i : i + subseqLen - 1) - mu(i));
    trend(i) = sum(z .* slope);
    invnorm(i) = 1 ./ norm(z - slope * trend(i));
end

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
df = [0; (1/2)*(timeSeries(1 + subseqLen : end) - timeSeries(1 : end - subseqLen))];
dg = [0; (timeSeries(1 + subseqLen : end) - mu(2 : end)) + (timeSeries(1 : n - subseqLen) - mu(1 : end - 1))];

matrixProfileIdx = NaN(length(mu), 1);

% The terms row and diagonal here refer to a hankel matrix representation of a time series
% This uses scaled cross correlation (scaled to allow the use of movstd) as an intermediate quantity for performance reasons.
% It is later reduced to z-normalized euclidean distance.
for diag = minlag + 1 : n - subseqLen + 1
    cov_ = sum((timeSeries(diag : diag + subseqLen - 1) - mu(diag)) .* (timeSeries(1 : subseqLen) - mu(1)));
    for row = 1 : n - subseqLen - diag + 2
        cov_ = cov_ + df(row) * dg(row + diag - 1) + df(row + diag - 1) * dg(row);
        corr_ = (cov_  - trend(row) * trend(row + diag - 1)) * invnorm(row) * invnorm(row + diag - 1);
        if corr_ > matrixProfile(row)
            matrixProfile(row) = corr_;
            matrixProfileIdx(row) = row + diag - 1;
        end
        if corr_ > matrixProfile(row + diag - 1)
            matrixProfile(row + diag - 1) = corr_;
            matrixProfileIdx(row + diag - 1) = row;
        end
    end
end


% This is the correct ordering. findMotifsDiscords uses a correlation based
% updates to avoid it being problematically slow if called in a loop.
[motifsIdx, mpAugmented] = findMotifsDetrended(timeSeries, mu, invnorm, matrixProfile, matrixProfileIdx, trend, slope, subseqLen, 3, 10, minlag, 2);
[discordsIdx] = findDiscords(mpAugmented, 3, minlag);
timeSeries(nanIdx) = NaN;
matrixProfile = (1 - min(1, matrixProfile, 'includenan'));

% Since we are plotting without explicit detrending, the scale factor has to be
% computed for the 'normal' case.
% invn_ = fastinvn(timeSeries, mu, subseqLen);
% mpgui.launchGui(timeSeries, mu, invn_, matrixProfile, motifsIdx, discordsIdx, subseqLen);


if transposed_  % matches the profile and profile index but not the motif or discord index to the input format
    matrixProfile = transpose(matrixProfile);
    matrixProfileIdx = transpose(matrixProfileIdx);
end

end

function [discordIdx] = findDiscords(matrixProfile, discordCount, exclusionLen)
% This function assumes a correlation based implementation of matrixProfile
% and that any excluded regions have been set to NaN, nothing else.
discordIdx = zeros(discordCount, 1);
maxCount = min(2 * discordCount * exclusionLen, numel(~isnan(matrixProfile)));
[corr_, idx] = mink(matrixProfile, maxCount);
for i = 1 : discordCount
    firstValid = idx(find(~isnan(corr_), 1));
    if isempty(firstValid)
        discordIdx = discordIdx(1 : i - 1);
        break;
    end
    discordIdx(i) = firstValid;
    % this assumes the chosen index is part of the exclusion length, so we
    % have up to exclusionLen - 1 excluded indices on each side. If exclusionLen is
    % greater than or equal to subsequence/motif length, then discords
    % cannot partially overlap.
    corr_(abs(idx - discordIdx(i)) < exclusionLen) = NaN;
end
end

function [motifIdxs, matrixProfile] = findMotifsDetrended(timeSeries, mu, invnorm, matrixProfile, profileIndex, trendComp, trend, subseqLen, motifCount, neighborCount, exclusionLen, radius)
% This is adapted match the output of some inline code written by Michael Yeh
% to find the top k motifs in a time series. Due to some bug fixes I applied, the two
% may return slightly different results, although they almost always agree on the top 2 motif pairs. 


motifIdxs = NaN(neighborCount, motifCount);
padLen = 2^nextpow2(length(timeSeries));
crosscov = @(idx) ifft(fft(timeSeries, padLen) .* conj(fft((timeSeries(idx : idx + subseqLen - 1) - mu(idx) - trendComp(idx) * trend) .* invnorm(idx), padLen)), 'symmetric');


for i = 1 : motifCount
    [corr_, motIdx] = max(matrixProfile);
    if ~isfinite(corr_)
        break;
    end
    % order subsequence motif pair as [time series index of 1st appearance, time series index of 2nd appearance]
    motifIdxs(1 : 2, i) = [min(motIdx, profileIndex(motIdx)), max(motIdx, profileIndex(motIdx))];
    [corrProfile] = crosscov(motIdx);
    corrProfile = min(1, corrProfile(1 : length(timeSeries) - subseqLen + 1) .* invnorm, 'includenan');
    % This uses correlation instead of normalized Euclidean distance,
    % because it's easier to work with and involves fewer
    % operations.
    
    corrProfile(isnan(matrixProfile)) = NaN;
    if exclusionLen > 0
        for j = 1 : 2
            exclRangeBegin = max(1, motifIdxs(j, i) - exclusionLen + 1); 
            exclRangeEnd = min(length(matrixProfile), motifIdxs(j, i) + exclusionLen - 1);
            corrProfile(exclRangeBegin : exclRangeEnd) = NaN;
        end
    end
    
    for j = 3 : neighborCount + 2
        [neighborCorr, neighbor] = max(corrProfile);
        if  ~isfinite(neighborCorr) || neighborCorr <= 0 || (1 - neighborCorr >= radius * (1 - corr_)) 
            break;
        end
        motifIdxs(j, i) = neighbor;
        if exclusionLen > 0
            exclRangeBegin = max(1, neighbor - exclusionLen + 1); 
            exclRangeEnd = min(length(matrixProfile), neighbor + exclusionLen - 1);
            corrProfile(exclRangeBegin : exclRangeEnd) = NaN;
            corrProfile(exclRangeBegin : exclRangeEnd) = NaN;
        end
    end
    % Matlab allows NaNs to be ignored from min/max calculations, so this
    % is a Matlab specific way to iteratively add excluded ranges of elements. 
    matrixProfile(isnan(corrProfile)) = NaN;
end
end

function [invn] = fastinvn(ts, mu, sublen)
% This is a simple variation on Welford's method. 
% This version still results in some minor cancellation and could be
% improved. It isn't prone to anything disastrous. 
invn = zeros(length(mu),1);
invn(1) = sum((ts(1:sublen)-mu(1)).^2);
for i = 2:length(invn)
    invn(i) = invn(i-1) + ((ts(i-1) - mu(i-1)) + (ts(i+sublen-1) - mu(i))) * (ts(i+sublen-1)-ts(i-1)); 
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
