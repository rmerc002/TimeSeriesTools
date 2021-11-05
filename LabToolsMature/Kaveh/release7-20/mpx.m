function [matrixProfile, matrixProfileIdx, mu, invnorm, preempted] = mpx(timeSeries, minlag, subseqLen, updateHandle)

% Code and skew symmetric update formulas by Kaveh Kamgar. They originated
% as a modification to earlier work by Yan Zhu.
% This version outputs a profile based on maximal correlation. It is mostly
% used in conjunction with mpx_anytime. 
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

preempted = false;
preemptible = false;

if nargin == 4
    preemptible = true;
elseif nargin ~= 3
    error('incorrect number of input arguments');
elseif ~isvector(timeSeries)
    error('first argument must be a 1D vector');
elseif ~(isfinite(subseqLen) && floor(subseqLen) == subseqLen) || (subseqLen < 2)
    error('subsequence length must be an integer value between 2 and the length of the timeSeries');
elseif minlag <= 0 || round(minlag) ~= minlag 
    error('Minimum lag/offset must be a positive integer');
elseif length(timeSeries) - minlag < subseqLen
    error('subsequence length + minimum lag/offset is too large');
end

transposed_ = isrow(timeSeries);
if transposed_
    timeSeries = transpose(timeSeries);
end

nanSubseqs = find(~isfinite(movsum(timeSeries, [0 subseqLen-1], 'Endpoints', 'discard')));
matrixProfile = repmat(-1, n - subseqLen + 1, 1);
matrixProfile(nanSubseqs) = NaN;
timeSeries(isnan(timeSeries)) = 0;
mu = moving_mean(timeSeries, subseqLen);
invnorm = fastinvn(timeSeries, mu, subseqLen);
invnorm(nanSubseqs) = NaN;

% The difference equations df and dg are slightly difficult to derive directly, but they arise naturally taking the difference
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

matrixProfileIdx = NaN(n - subseqLen + 1, 1);

% The terms row and diagonal here refer to a hankel matrix representation of a time series
% This uses scaled partial cross correlation as an intermediate quantity for performance reasons.
% It is later reduced to z-normalized euclidean distance.
for diag = minlag + 1 : n - subseqLen + 1
    drawnow;
    if preemptible
        sh = updateHandle.shouldHalt();
        if sh
            fprintf('abandoning current iteration (time series length: %d subsequence length: %d)\n', n, subseqLen);
            preempted = true;
            break;
        end
    end
    cov_ = (sum((timeSeries(diag : diag + subseqLen - 1) - mu(diag)) .* (timeSeries(1 : subseqLen) - mu(1))));
    for row = 1 : n - subseqLen - diag + 2
        cov_ = cov_ + df(row) * dg(row + diag - 1) + df(row + diag - 1) * dg(row);
        corr_ = cov_ * invnorm(row) * invnorm(row + diag - 1);
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


if transposed_  % matches the profile and profile index but not the motif or discord index to the input format
    matrixProfile = transpose(matrixProfile);
    matrixProfileIdx = transpose(matrixProfileIdx);
end

end
