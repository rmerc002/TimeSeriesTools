function [mpa, mpb, mpia, mpib] = abbaJoindetrend(a, b, subseqLen, selfJoin)
% References
% Yan Zhu, et al, Matrix Profile II: Exploiting a Novel Algorithm and GPUs to break the one Hundred Million Barrier for Time Series Motifs and Join
% Zachary Zimmerman, et al, Scaling Time Series Motif Discovery with GPUs: Breaking the Quintillion Pairwise Comparisons a Day Barrier. (pending review)
% Philippe Pebay, et al, Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments
% Takeshi Ogita, et al, Accurate Sum and Dot Product
%

% difference equations have 0 as their first entry here to simplify index
% calculations slightly. Alternatively, it's also possible to swap this to the last element
% and reorder the comparison step (or omit on the last step). This is a
% special case when comparing a single time series to itself. The more general
% case with time series A,B can be computed using difference equations for
% each time series.

if nargin ~= 4
    error('incorrect number of input arguments');
elseif ~isvector(a) || ~isvector(b)
    error('first argument must be a 1D vector');
elseif ~(isfinite(subseqLen) && floor(subseqLen) == subseqLen) || (subseqLen < 2) || (subseqLen > length(a)) || (subseqLen > length(b))
    error('subsequence length must be an integer value between 2 and the length of the timeSeries');
end

transposed_ = isrow(a);
if transposed_
    if iscol(b)
        error('This function does not accept different input formats for a and b, because it leads to ambiguous output formatting (and typically downstream errors)');
    end
    a = transpose(a);
    b = transpose(b);
end

slope = fittingLine(subseqLen);

[a, mua, invnorma, dfa, dga, trenda] = preprocess(a, subseqLen, slope);
[b, mub, invnormb, dfb, dgb, trendb] = preprocess(b, subseqLen, slope);

mpa = repmat(-1, length(a) - subseqLen + 1, 1);
mpb = repmat(-1, length(b) - subseqLen + 1, 1);
mpa(~isfinite(invnorma)) = nan;
mpb(~isfinite(invnormb)) = nan;
mpia = NaN(length(a) - subseqLen + 1, 1);
mpib = NaN(length(b) - subseqLen + 1, 1);

% I'm not suggesting this is always robust. NaNs compare false against
% everything when using ordered comparisons. 
% This means 
% (nan < val)  is false
% (nan > val)  is false
% (nan == val) is false      

% On min/max ops, it is much less consistent. It's a matter of how they're
% implemented.
% Most implementations favor implementing min as 
% if a < b return a else return b
% and max as
% if a > b return a else return b

% but a few of them (notably CUDA and Matlab's default environment) differ
% For matlab, the default is omitnan. I'm sticking with update on true
% comparison, but don't assume this is always robust. It is very subtle,
% and it can harbor subtle bugs. I'm doing it because it was requested. 

[mpa, mpb, mpia, mpib] = symmABJoin(a, b, mua, mub, invnorma, invnormb, dfa, dga, dfb, dgb, trenda, trendb, mpa, mpb, mpia, mpib, subseqLen, selfJoin);
[mpb, mpa, mpib, mpia] = symmABJoin(b, a, mub, mua, invnormb, invnorma, dfb, dgb, dfa, dga, trendb, trenda, mpb, mpa, mpib, mpia, subseqLen, selfJoin);

mpa = max(0, 1 - mpa, 'includenan');
mpb = max(0, 1 - mpb, 'includenan');

if transposed_
    mpa = transpose(mpa);
    mpb = transpose(mpb);
    mpia = transpose(mpia);
    mpib = transpose(mpib);
end

end

function [mpa, mpb, mpia, mpib] = symmABJoin(a, b, mua, mub, invnorma, invnormb, dfa, dga, dfb, dgb, trendCompa, trendCompb, mpa, mpb, mpia, mpib, w, selfJoin)
subseqCounta = length(a) - w + 1;
subseqCountb = length(b) - w + 1;
for ia = 1 : subseqCounta
    mx = min(subseqCounta - ia + 1, subseqCountb);
    cov_ = sum((a(ia : ia + w - 1) - mua(ia)) .* (b(1 : w) - mub(1)));
    for ib = 1 : mx
        if selfJoin && ia-1 < w/2
            continue
        end
        
        cov_ = cov_ + dfa(ib + ia - 1) * dgb(ib) + dga(ib + ia - 1) * dfb(ib);
        corr_ = (cov_ - trendCompa(ib + ia - 1) * trendCompb(ib)) * invnorma(ib + ia - 1) * invnormb(ib);
        if corr_ > mpa(ib + ia - 1)
            mpa(ib + ia - 1) = corr_;
            mpia(ib + ia - 1) = ib;
        end
        if corr_ > mpb(ib)
            mpb(ib) = corr_;
            mpib(ib) = ib + ia - 1; 
        end
    end
end
end


function slope = fittingLine(subseqLen)
center = ceil(subseqLen / 2);
% These are just variations of the classic formula for the sum of squares from 1 to center - 1
if mod(subseqLen,2) == 0
   slope = (-center + 0.5 : center - 0.5)' ./ sqrt((center * (center + 1) * (2 * center + 1)) / 3 - center * (center + 1) + 0.5 * center);
else
   slope = (-center + 1 : center - 1)' ./ sqrt((center - 1) * center * (2 * center - 1) / 3);
end
end

% factored this part out because I was tired of duplicates
function [timeSeries, mu, invnorm, df, dg, trend] = preprocess(timeSeries, subseqLen, slope)

nanSubseq = ~isfinite(movsum(timeSeries, [0 subseqLen-1], 'Endpoints', 'discard'));
timeSeries(~isfinite(timeSeries)) = 0;
mu = moving_mean(timeSeries, subseqLen);
subseqCount = length(mu);
df = [0; (1/2)*(timeSeries(1 + subseqLen : end) - timeSeries(1 : end - subseqLen))];
dg = [0; (timeSeries(1 + subseqLen : end) - mu(2 : end)) + (timeSeries(1 : end - subseqLen) - mu(1 : end - 1))];

trend = zeros(subseqCount, 1);
invnorm = zeros(subseqCount, 1);

for i = 1 : subseqCount
    z = (timeSeries(i : i + subseqLen - 1) - mu(i));
    trend(i) = sum(z .* slope);
    if ~nanSubseq(i)
        invnorm(i) = 1 ./ norm(z - slope * trend(i));
    else
        invnorm(i) = nan;
    end
end

end


function [ res ] = moving_mean(a,w)
% moving mean over sequence a with window length w
% based on Ogita et. al, Accurate Sum and Dot Product

% This moving mean function is used to reduce accumulated rounding error (in spite of the fact that motif
% discovery isn't that sensitive to them).
% Using a typical summation algorithm to compute the mean for each window independently, the
% resulting errors attributable to accumulated roundoff form a not.0
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

