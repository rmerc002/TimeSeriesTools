function [matrixProfile, matrixProfileIdx, motifsIdx, discordsIdx] = mpx_v2(timeseries, minlag, subseqlen, plot_output)
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

% Code and update formulas are by Kaveh Kamgar. 
% GUI and top k motif critera are based on but not identical to some code
% by Michael Yeh. Implementation details are specified above.
%
% Additional References
% Yan Zhu, et al, Matrix Profile II: Exploiting a Novel Algorithm and GPUs to break the one Hundred Million Barrier for Time Series Motifs and Join
% Zachary Zimmerman, et al, Scaling Time Series Motif Discovery with GPUs: Breaking the Quintillion Pairwise Comparisons a Day Barrier. (pending review)
% Philippe Pebay, et al, Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments
% Takeshi Ogita, et al, Accurate Sum and Dot Product

subcount = length(timeseries) - subseqlen + 1;

% difference equations have 0 as their first entry here to simplify index
% calculations slightly. Alternatively, it's also possible to swap this to the last element
% and reorder the comparison step (or omit on the last step). This is a
% special case when comparing a single time series to itself. The more general
% case with time series A,B can be computed using difference equations for
% each time series.

if nargin == 3
    plot_output = true;
elseif nargin ~= 4
    error('incorrect number of input arguments');
elseif ~isvector(timeseries)
    error('first argument must be a 1D vector');
elseif ~(isfinite(subseqlen) && floor(subseqlen) == subseqlen) || (subseqlen < 2) || (subcount < 2) 
    error('subsequence length must be an integer value between 2 and the length of the timeseries');
end

transposed_ = isrow(timeseries);
if transposed_
    timeseries = transpose(timeseries);
end

nanmap = find(~isfinite(movsum(timeseries, [0 subseqlen-1], 'Endpoints', 'discard')));
nanIdx = find(isnan(timeseries));
timeseries(nanIdx) = 0;
% We need to remove any NaN or inf values before computing the moving mean,
% because it uses an accumulation based method. We add additional handling
% elsewhere as needed.


mu = moving_mean(timeseries, subseqlen);
mus = moving_mean(timeseries, subseqlen - 1);
invnorm = zeros(subcount, 1);

% Methods fail more often here, including the builtin matlab movstd.
% This uses a simple approach, because it's better at detecting constant
% regions, which should be completely ignored.
for i = 1 : subcount
    invnorm(i) = 1 ./ norm(timeseries(i : i + subseqlen - 1) - mu(i));
end

invnorm(nanmap) = NaN;
invnorm(~isfinite(invnorm)) = NaN;

% This method requires 4 arrays for rank 1 updates. The previous method used
% 2 arrays in the case of iterating over a single time series. This seems
% to have better stability when dealing with time series with missing data.
% It still loses a number of digits of precision in very badly conditioned
% data.

% References
% Phillippe Pebay, Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments
% Knuth, semi-numerical algorithms, vol 2

dr_bwd = [0; timeseries(1 : subcount - 1) - mu(1 : subcount - 1)];
dc_bwd = [0; timeseries(1 : subcount - 1) - mus(2 : subcount)];
dr_fwd = timeseries(subseqlen : end) - mu(1 : subcount);
dc_fwd = timeseries(subseqlen : end) - mus(1 : subcount);
matrixProfile = repmat(-1, subcount, 1);
matrixProfile(nanmap) = NaN;
matrixProfileIdx = NaN(subcount, 1);

for diag = minlag + 1 : subcount
    cov_ = sum((timeseries(diag : diag + subseqlen - 1) - mu(diag)) .* (timeseries(1 : subseqlen) - mu(1)));
    for row = 1 : subcount - diag + 1
        col = diag + row - 1;
        if row > 1
            cov_ = cov_ - dr_bwd(row) * dc_bwd(col) + dr_fwd(row) * dc_fwd(col);
        end
        corr_ = cov_ * invnorm(row) * invnorm(col);
        if corr_ > matrixProfile(row)
            matrixProfile(row) = corr_;
            matrixProfileIdx(row) = col;
        end
        if corr_ > matrixProfile(col)
            matrixProfile(col) = corr_;
            matrixProfileIdx(col) = row;
        end
    end
end


% updated to pick outliers independent of motifs
% this means they might overlap or whatever. Compute more of them if you
% need to in order to avoid this. The other methods introduce bias. 
[discordsIdx] = findDiscords(matrixProfile, minlag);

[motifsIdx] = findMotifs(timeseries, mu, invnorm, matrixProfile, matrixProfileIdx, subseqlen, minlag);


% Max caps anything that rounds inappropriately. This can hide potential
% issues at times, but it's fairly rare. 
matrixProfile = sqrt(max(0, 2 * subseqlen * (1 - matrixProfile), 'includenan'));
timeseries(nanIdx) = NaN;

if plot_output
    gui = mpgui_ed(timeseries, subseqlen);
    % first pair is plottted alongside data
    gui.plotProfile(matrixProfile);
    if isfinite(motifsIdx(1, 1)) && ~isnan(motifsIdx(2, 1))
        gui.plotData(motifsIdx(1, 1), motifsIdx(2, 1));
        motiftitle = sprintf('The first motif pair is located at %d and %d.', motifsIdx(1,1), motifsIdx(2,1));
        gui.plotMotif(gui.motifAx1, motifsIdx(1, 1), motifsIdx(2:end, 1), motiftitle);
        if isfinite(motifsIdx(1, 2))
            motiftitle = sprintf('The second motif pair is located at %d and %d.', motifsIdx(1,2), motifsIdx(2,2));
            gui.plotMotif(gui.motifAx2, motifsIdx(1, 2), motifsIdx(2:end, 2), motiftitle);
        end
        if isfinite(motifsIdx(1, 3))
            motiftitle = sprintf('The third motif pair is located at %d and %d.', motifsIdx(1,3), motifsIdx(2,3));
            gui.plotMotif(gui.motifAx3, motifsIdx(1, 3), motifsIdx(2:end, 3), motiftitle);
        end
        gui.plotDiscords(discordsIdx, sprintf(['The top three discords ', '%d(blue), %d(red), %d(green)'], discordsIdx(1), discordsIdx(2), discordsIdx(3)));
    else
        gui.plotData();
    end
    
    gui.drawgui;
    
end

if transposed_  % matches the profile and profile index but not the motif or discord index to the input format
    matrixProfile = transpose(matrixProfile);
    matrixProfileIdx = transpose(matrixProfileIdx);
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
