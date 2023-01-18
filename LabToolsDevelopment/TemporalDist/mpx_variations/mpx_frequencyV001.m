function [frequencyProfile] = mpx_frequencyV001(timeseries, minlag, subseqlen, numThresholds)
% matrixProfile - distance between each subsequence timeseries(i : i + subseqLen - 1)
%                 for i = 1 .... length(timeseries) - subseqLen + 1, with the condition
%                 that the nearest neighbor of the ith subsequence is at least minlag steps
%                 away from it.
%
% matrixProfileIdx - If matrixProfileIdx(i) == j, then |i - j| >= minlag
%                    and MatrixProfile(i) == norm(zscore(subseq(i), 1) - norm(zscore(subseq(j), 1));
%
% motifsIdx - A 10 x 3 matrix containing the indices of the top motifs and
%             their neighbors. The first two entries in each column
%             indicate a pair. The remaining entries in that column are their
%             neighbors. This implementation picks motifs from best to
%             worst. It picks neighbors using a greedy strategy with an
%             excluded region applied after each pick to all subsequent
%             choices of motifs and neighbors. Outliers are picked first,
%             followed by the top motif and its neighbors, followed by the
%             second and third motifs in the same manner. This avoids
%             picking outliers as motifs on short time series.
%
%             The distance between a motif and its neighbor must be less than
%             two times the distance between a motif and its nearest
%             neighbor. This is not strict enough to guarantee any causal
%             relation, and in cases of weak matches, neighbors may vary
%             independently from the first pair (that is they may be
%             uncorrelated).
%
% discordsIdx - A 3 x 1 colum indicating the furthest outliers by maximal
%               normalized euclidean distance. These are chosen prior to
%               motifs in the current implementation.
%
%

% Edited 2/7/20: Changed update formula again due to numerical cancellation
% issues in some extreme cases, particularly ill conditioned time series
% with missing data. Indexing feature is still not robust to very poorly
% conditioned data.
%
% Edited: 12/24/20: Added missing data features. These are still somewhat
% experimental. They check for any non-normalizable windows, mark them, 
% and apply noise to any sequences of constants. The extra check of update 
% or no update is implicit, since the comparison is false if either operand
% is nan. This is not advisable in all environments.
%
% The difference equations formula was also adjusted to omit the leading 
% zero, because this casues problems in tiled implementations when starting
% from some point other than the beginning of the time series.

% Code and update formulas are by Kaveh Kamgar.
% GUI and top k motif critera are based on but not identical to some code
% by Michael Yeh. Implementation details are specified above.
%
% Additional References
% Yan Zhu, et al, Matrix Profile II: Exploiting a Novel Algorithm and GPUs to break the one Hundred Million Barrier for Time Series Motifs and Join
% Zachary Zimmerman, et al, Scaling Time Series Motif Discovery with GPUs: Breaking the Quintillion Pairwise Comparisons a Day Barrier. (pending review)
% Philippe Pebay, et al, Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments
% Takeshi Ogita, et al, Accurate Sum and Dot Product

initial_subcount = length(timeseries) - subseqlen + 1;

% difference equations have 0 as their first entry here to simplify index
% calculations slightly. Alternatively, it's also possible to swap this to the last element
% and reorder the comparison step (or omit on the last step). This is a
% special case when comparing a single time series to itself. The more general
% case with time series A,B can be computed using difference equations for
% each time series.

% if nargin == 3
%     plot_output = true;
% elseif nargin ~= 4
%     error('incorrect number of input arguments');
% elseif ~isvector(timeseries)
%     error('first argument must be a 1D vector');
% elseif ~(isfinite(subseqlen) && floor(subseqlen) == subseqlen) || (subseqlen < 2) || (initial_subcount < 2)
%     error('subsequence length must be an integer value between 2 and the length of the timeseries');
% end

transposed_ = isrow(timeseries);
if transposed_
    timeseries = transpose(timeseries);
end

[ts, isvalidwindow, first, last] = find_valid_windows(timeseries, subseqlen);
subcount = length(ts) - subseqlen + 1;
if subcount < minlag
    warning("no valid comparisons could be performed");
    frequencyProfile = inf(numThresholds, subcount);
    motifsIdx = [];
    discordsIdx = [];
    return;
end

mu = moving_mean(ts, subseqlen);
mus = moving_mean(ts(2:end-1), subseqlen - 1);
invnorm = zeros(subcount, 1);


% Methods fail more often here, including the builtin matlab movstd.
% This uses a simple approach, because it's better at detecting constant
% regions, which should be completely ignored.
for i = 1 : subcount
    invnorm(i) = 1 ./ norm(ts(i : i + subseqlen - 1) - mu(i));
end

if any(~isfinite(invnorm) & isvalidwindow)
    % This is an indicator that all of the preconditioners failed. It may
    % or may not impact the result.
    warning('non finite values encountered when computing reciprocal per window norm.');
end

invnorm(~isvalidwindow) = NaN;
invnorm(~isfinite(invnorm)) = NaN;

% This method requires 4 arrays for rank 1 updates. The previous method used
% 2 arrays in the case of iterating over a single time series. This seems
% to have better stability when dealing with time series with missing data.
% It still loses a number of digits of precision in very badly conditioned
% data.

% References
% Phillippe Pebay, Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments
% Knuth, semi-numerical algorithms, vol 2

% This was changed to remove the leading zero, as the leading zero
% makes for a bad reference. It doesn't work if you need to tile the
% problem and start from row > 1 rather than row == 1.
dr_bwd = ts(1 : subcount - 1) - mu(1 : subcount - 1);
dc_bwd = ts(1 : subcount - 1) - mus;
dr_fwd = ts(subseqlen + 1 : end) - mu(2 : subcount);
dc_fwd = ts(subseqlen + 1 : end) - mus;
frequencyProfile = zeros(numThresholds, subcount);
% matrixProfile(~isvalidwindow) = 0;
matrixProfileIdx = zeros(subcount,1);

thresholdVals = linspace(0,1,numThresholds);

for diag = minlag + 1 : subcount
    cov_ = sum((ts(diag : diag + subseqlen - 1) - mu(diag)) .* (ts(1 : subseqlen) - mu(1)));
    for row = 1 : subcount - diag + 1
        col = diag + row - 1;
        if row > 1
            cov_ = cov_ - dr_bwd(row-1) * dc_bwd(col-1) + dr_fwd(row-1) * dc_fwd(col-1);
        end
        corr_ = cov_ * invnorm(row) * invnorm(col);

        threshIndices = corr_ > thresholdVals;

        frequencyProfile(threshIndices, row) = frequencyProfile(threshIndices, row) + 1;
        frequencyProfile(threshIndices, col) = frequencyProfile(threshIndices, col) + 1;

    end
end

for ii = 1:numThresholds
    frequencyProfile(ii,:) = frequencyProfile(ii,:)./length(frequencyProfile(ii,:));
end

% for ii = 1:numThresholds
%     frequencyProfile(ii,:) = frequencyProfile(ii,:)./max(frequencyProfile(ii,:));
% end

figure;
imagesc(frequencyProfile);

set(gca,'YDir','normal');
end

function [timeseries, isvalidwindow, first, last] = find_valid_windows(timeseries, subseqlen)

% We're going to check for non-normalizable windows here.

% Check for windows containing non-finite data.
subcount = length(timeseries) - subseqlen + 1;
isfinitedata = isfinite(timeseries);
isvalidwindow = movsum(isfinitedata, [0 subseqlen-1], 'Endpoints', 'discard') == subseqlen;

% Typically these come in streams, so we zero them rather than use the
% windowed mean

timeseries(~isfinitedata) = 0;

% Now we need to check for sequences of constants, since these represent
% singularities, which greatly impact stability

issingularity = false(subcount, 1);
% find points with 2 identical consecutive points
for i = 1 : subcount
    % it doesn't matter whether any of these are attributable
    % to zeroing, as this is just a conditioning issue
    issingularity(i) = all(timeseries(i + 1 : i + subseqlen - 1) == timeseries(i));
end

isvalidwindow(issingularity) = false;

% Trim cases where we have bad leading or trailing windows
first = find(isvalidwindow, 1);
last = find(isvalidwindow, 1, 'last');
timeseries = timeseries(first : last + subseqlen - 1);
isvalidwindow = isvalidwindow(first : last);

if isempty(timeseries) || last == first
    return;
end

issingularity = issingularity(first : last);
singularities = find(issingularity);

i = 1;
while i < length(singularities)
    j = i;
    while j < length(singularities)
        % finding the end
        if singularities(j + 1) ~= singularities(j) + 1
            break;
        end
        j = j + 1;
    end
    % find an appropriate scale factor based on the constant occupying
    % this region.
    v = abs(timeseries(singularities(i)));
    if v == 0 
        %just use a standard deviation of 1.
        % It might seem better to look at surrounding windows, but at 
        % that point, we have to recursively check those windows.
        scale = 1;
    elseif v < 1
        % 0 < v < 1, so natural log is negative
        p = 1 + log(v);
        if p == 0
            % 1/2^52 on most systems
            scale = eps(1);
        else
            scale = 1 / abs(p);
        end
    else
        scale = 1 + log(v);
    end
   
    % 0.01 is just to ensure this doesn't perturb leading digits
    % in case this data provides a partial window for some normalizable
    % section. Feel free to make it smaller or larger. 
    % The data variance dependent portion is handled by "scale".
    timeseries(i : j + subseqlen - 1) = 0.01 * scale * randn(j + subseqlen - i, 1);
    
    % If any kind of error happens here, this is the point where we give
    % up, as it becomes very difficult to achieve predictable results.
    if any(~isfinite(timeseries(i : j + subseqlen - 1)))
        error(fprintf('unable to fix missing spanning %d to %d\n', i, j));
    end
    
    i = j + 1;
    
end

end


function [discordIdx] = findDiscords(matrixProfile, exclusionLen)
% This function assumes a correlation based implementation of matrixProfile
% and that any excluded regions have been set to NaN, nothing else.

% Changed on 1/29/20
% Making motif discovery and discord discovery co-dependent introduces odd
% issues. I am avoiding any sharing of excluded regions at this point,
% since it's not clear which should come first. Computing 3 or k or
% whatever motifs prior to discords was definitely not correct, but I'm not
% sure this is a better idea.

discordCount = 3;
discordIdx = zeros(discordCount, 1);
[~, idx] = sort(matrixProfile);

for i = 1 : discordCount
    f = find(isfinite(matrixProfile(idx)) & (matrixProfile(idx) > -1), 1);
    if isempty(f) || (matrixProfile(f) == -1)
        discordIdx(i : end) = NaN;
        break;
    end
    discordIdx(i) = idx(f);
    exclRangeBegin = max(1, discordIdx(i) - exclusionLen + 1);
    exclRangeEnd = min(length(matrixProfile), discordIdx(i) + exclusionLen - 1);
    matrixProfile(exclRangeBegin : exclRangeEnd) = NaN;
end
end


function [motifIdxs, matrixProfile] = findMotifs(timeseries, mu, invnorm, matrixProfile, profileIndex, subseqLen, exclusionLen)
% This is adapted match the output of some inline code written by Michael Yeh
% to find the top k motifs in a time series. Due to some bug fixes I applied, the two
% may return slightly different results, although they almost always agree on the top 2 motif pairs.


% +2 allows us to pack the original motif pair with its neighbors.
% Annoyingly if it's extended, matlab will pad with zeros, leading to some
% rather interesting issues later.
motifCount = 3;
radius = 2;
neighborCount = 10;

motifIdxs = NaN(neighborCount + 2, motifCount);

% removing use of ffts since this formula is easier to analyze and conv2
% is pretty well optimized
crosscov = @(idx) conv2(timeseries, (timeseries(idx + subseqLen - 1 : -1 : idx) - mu(idx)) .* invnorm(idx), 'valid');

for i = 1 : motifCount
    [corr_, motIdx] = max(matrixProfile);
    % -1 means this is maximally negatively correlated and was therefore
    % never updated
    if ~isfinite(corr_) || corr_ == -1
        break;
    end
    % order subsequence motif pair as [time series index of 1st appearance, time series index of 2nd appearance]
    motifIdxs(1 : 2, i) = [min(motIdx, profileIndex(motIdx)), max(motIdx, profileIndex(motIdx))];
    [corrProfile] = crosscov(motIdx);
    corrProfile = min(1, corrProfile(1 : length(timeseries) - subseqLen + 1) .* invnorm, 'includenan');
    % This uses correlation instead of normalized Euclidean distance, because it's easier to work with and involves fewer operations.
    
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
        % If you want to put a reasonable global bound on it,
        % set the if statement to also skip anything where neighborCorr <= 0 as well.
        % This was reverted to match earlier code, which did not enforce such a restriction.
        %
        if  ~isfinite(neighborCorr) || ((1 - neighborCorr) >= radius * (1 - corr_))
            break;
        end
        motifIdxs(j, i) = neighbor;
        if exclusionLen > 0
            exclRangeBegin = max(1, neighbor - exclusionLen + 1);
            exclRangeEnd = min(length(matrixProfile), neighbor + exclusionLen - 1);
            corrProfile(exclRangeBegin : exclRangeEnd) = NaN;
        end
    end
    % Matlab allows NaNs to be ignored from min/max calculations by default, so this
    % is a Matlab specific way to iteratively add excluded ranges of elements.
    matrixProfile(isnan(corrProfile)) = NaN;
end
end


function [ res ] = moving_mean(a,w)
% moving mean over sequence a with window length w
% based on Ogita et. al, Accurate Sum and Dot Product

% A major source of rounding error is accumulated error in the mean values, so we use this to compensate.
% While the error bound is still a function of the conditioning of a very long dot product, we have observed
% a reduction of 3 - 4 digits lost to numerical roundoff when compared to older solutions.

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
