function [plato, platoIndex, platoContrast] = RelativeFrequencyContrastProfile(positiveTS, negativeTS, m, maxFreq, forcePlot)
    %%% Input:
    %%%   positiveTS: A time series containing at least two instances of a desired behavior
    %%%   negativeTS: A time series containing zero instances of a desired behavior
    %%%   m: The approximate window size of the desired behavior
    %%%   maxFreq: Number of nearest neighbors to search. Does not affect
    %%%     time complexity
    %%%   forcePlot: (optional) Display the following two plots
    %%%
    %%% Output:
    %%%   plato: The subsequence that most distinguishes
    %%%          positiveTS from negativeTS. Determined based on
    %%%          maximum of root mean squared contrast for each time index.
    %%%   RFCP: Contrast Profile, which indicates subsequence within
    %%%       positiveTS that are less conserved in negativeTS
    %%%   RFMP_AA_Indices: K nearest neighbor indices within positiveTS

    initial__pos_subcount = length(positiveTS) - m + 1;
    initial__neg_subcount = length(negativeTS) - m + 1;

    if nargin == 4
        forcePlot = false;
    elseif nargin ~= 5
        error('incorrect number of input arguments');
    elseif ~isvector(positiveTS) || ~isvector(negativeTS)
        error('first and second arguments must be a 1D vector');
    elseif ~(isfinite(m) && floor(m) == m) || (m < 2) || (initial__pos_subcount < 2) || (initial__neg_subcount < 2)
        error('subsequence length must be an integer value between 2 and the length of the timeseries');
    end

%     tsPos(isnan(tsPos)) = nanmean(tsPos);
%
%     tsNeg(isnan(tsNeg)) = nanmean(tsNeg);

    %Change to row vector for internal use
    if size(positiveTS,1) == 1
        %save the data orientation for matching output format to input
        %give orientation priority to positiveTS
        positiveTS = positiveTS';
    end
    if size(negativeTS,1) == 1
        %give orientation priority to positiveTS, do not save negative
        %orientation if different than positive
        negativeTS = negativeTS';
    end



    %%% Matrix profile self-join using positiveTS
    [RFMP_AA, RFMP_AA_Indices] = RelativeFrequencyMatrixProfile(positiveTS, positiveTS, m, maxFreq);
    RFMP_AA = real(RFMP_AA);
    %TODO: need to figure out the bound of MASS_V2 in mpxAAKNN

    %%% Euclidean distance values above sqrt(2*m) are equivalent to
    %%%   anti-correlated values
    RFMP_AA = clipMatrixProfileAmplitude(RFMP_AA, m);
    %%% pad with NaN to make future comparisons between matrix profiles
%     padLength = length(positiveTS) - length(MP_AA) + 1;
%     MP_AA = [MP_AA;NaN(padLength,1)];

    %%% Matrix profile AB-join between positiveTS and negativeTS
    [RFMP_AB, RFMP_AB_Indices] = RelativeFrequencyMatrixProfile(positiveTS, negativeTS, m, maxFreq);
    RFMP_AB = real(RFMP_AB);
    %TODO: need to figure out the bound of MASS_V2 in mpxAAKNN

    %%% Euclidean distance values above sqrt(2*m) are equivalent to
    %%%   anti-correlated values
    RFMP_AB = clipMatrixProfileAmplitude(RFMP_AB, m);
    %%% pad with NaN to make future comparisons between matrix profiles
%     padLength = length(positiveTS) - length(MP_AB) + 1;
%     MP_AB = [MP_AB;NaN(padLength,1)];

    %%%If a non-overlapping candidates is less than maxFreq
    effectiveMaxFreq = min(size(RFMP_AA,1), size(RFMP_AB,1));
    if maxFreq > effectiveMaxFreq
        fprintf("Warning!: Choice of maxFreq=%d was too large for dataset, setting to %d\n",maxFreq, effectiveMaxFreq);
        RFMP_AA = RFMP_AA(1:effectiveMaxFreq,:);
        RFMP_AB = RFMP_AB(1:effectiveMaxFreq,:);
    end



    %%% Contrast Profile
    RFCP = RFMP_AB - RFMP_AA;
    %%% Normalize values to the range [0,1]
    RFCP = normalizeContrastProfileAmplitude(RFCP, m);

    %%% For Debugging
%     figure;
%     tiledlayout(3,1);
%     ax1 = nexttile();
%     plot(positiveTS);
%
%     ax2 = nexttile();
%     hold on;
%     plot(MP_AA(1,:));
%     plot(MP_AB(1,:));
%     hold off;
%
%     ax3 = nexttile();
%     hold on;
%     plot(CP(1,:));
%     plot(CP(7,:));
%     hold off;
%     linkaxes([ax1, ax2, ax3], 'x');

    %%% plato is the subsequence in positiveTS corresponding to index with
    %%%   largest contrast profile value
%     [maxContrastValues, platoIndices] = max(RFCP,[],2);
%     [maxContrastValue, KBest] = max(maxContrastValues);
%     platoIndex = platoIndices(KBest);
%     plato = tsPos(platoIndex:platoIndex+m-1);

    %%%Alternative Plato calculation
    %%%Maximum root mean squared of contrast values
    CPmean = zeros(1,size(RFCP,2));
    for ti = 1:size(RFCP,2)-m+1
        subsequence = positiveTS(ti:ti+m-1);
        if sum(isnan(subsequence)) > 0
            continue
        end
       CPmean(ti) = norm(RFCP(:,ti))/sqrt(maxFreq); %% root mean squared
    end
    [platoContrast, platoIndex] = max(CPmean);
    platoIndex = platoIndex;
    startIndex = platoIndex;
    endIndex = startIndex + m - 1;
    plato = positiveTS(startIndex:endIndex);


    %%% plato_twin is the nearest neighbor of plato within positieTS
    %%%   It is not necessarily the second largest contrast profile value
    platoTwinIndex = RFMP_AA_Indices(1, platoIndex);
    if platoTwinIndex+ m-1 > length(positiveTS)
        platoTwinIndex = length(positiveTS)-m+1;
    end
%     plato_twin = tsPos(platoTwinIndex:platoTwinIndex + m - 1);



    %%%%%%%%%%%%%
    %%% PLOTS %%%
    %%%%%%%%%%%%%
    if forcePlot == true
        visualizeMMPAB(positiveTS, 1-RFCP, RFMP_AA_Indices, [], [], [], m, 1:maxFreq, "KNN","");

        redColor = [0.73,0.05,0];
        greenColor = [0,0.73,0.41];
        blueColor = [0,0.29,0.73];
        grayColor = [0.75,0.75,0.75];
        lightBlueColor = [0.01, 0.83,0.99];
        platoColor = [129/255, 51/255, 144/255];
        platoTwinColor = [115/255, 170/255, 43/255];

        %%%%%%%%%%%%%%%%%%
        %%% PLATO plot %%%
        %%%%%%%%%%%%%%%%%%
        fig = figure('Name','Contrast Profile: Plato','NumberTitle','off');
        set(gcf, 'Position', [0,100,400,400]);
        annotation('textbox', [0, 0.75, 0, 0], 'string', 'Plato')
        annotation('textbox', [0, 0.35, 0, 0], 'string', 'Plato Twin')

        subsequenceIndices = [platoIndex, platoTwinIndex];
        colors = [platoColor; platoTwinColor];
        plotIndex = 0;
        inset = 0.9;
%         plot(0,0);
        hold on;
        for i = 1:length(subsequenceIndices)
            ti = subsequenceIndices(i);
            color = colors(i,:);
            tempTS = positiveTS(ti:ti+m-1);
            tempMin = min(tempTS);
            tempMax = max(tempTS);
            tempRange = max(1e-5, tempMax-tempMin);
            plot(-plotIndex+inset*(tempTS-tempMin)/tempRange,'Color',color);
            plotIndex = plotIndex + 1;
        end
        hold off;


        xlim([0,m]);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}Plato \\color{black}(top) and \\color[rgb]{%f,%f,%f}Plato Twin \\color{black}(bottom)",platoColor(1), platoColor(2), platoColor(3), platoTwinColor(1), platoTwinColor(2), platoTwinColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,m],'ytick',[], 'TickDir','out');
        box off;

        fprintf('plato index: %d, plato twin index: %d\n', platoIndex, platoTwinIndex);
    end
end

function [mp] = clipMatrixProfileAmplitude(mp,m)
    mp = min(sqrt(2*m),mp,'includenan'); %clip the anti-correlated values
    mp = real(mp);
end

function [cp] = normalizeContrastProfileAmplitude(cp, m)
    cp = cp/(sqrt(2*m)); %normalize between 0,1
    %discard negative values. These occur when a subsequence in T+ has a
    %closer nearest neighbor in T-. It's not a behavior of interest for
    %this research.
    cp = nanmax(0, cp);
end

function [matrixProfile, matrixProfileIdx, isvalidwindow, motifsIdx, discordsIdx] = mpxSelfJoin(timeseries, minlag, subseqlen, plot_output)
% (mpx_v3 titled from Kaveh Kamgar, 2020-12-24 12:58)
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

if nargin == 3
    plot_output = true;
elseif nargin ~= 4
    error('incorrect number of input arguments');
elseif ~isvector(timeseries)
    error('first argument must be a 1D vector');
elseif ~(isfinite(subseqlen) && floor(subseqlen) == subseqlen) || (subseqlen < 2) || (initial_subcount < 2)
    error('subsequence length must be an integer value between 2 and the length of the timeseries');
end

transposed_ = isrow(timeseries);
if transposed_
    timeseries = transpose(timeseries);
end

[ts, isvalidwindow, first, last] = find_valid_windows_SelfJoin(timeseries, subseqlen);
subcount = length(ts) - subseqlen + 1;
if subcount < minlag
    warning("no valid comparisons could be performed");
    matrixProfile = inf(subcount, 1);
    matrixProfileIdx = repmat(-1, subcount, 1);
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
matrixProfile = repmat(-1, subcount, 1);
matrixProfile(~isvalidwindow) = NaN;
matrixProfileIdx = NaN(subcount, 1);

for diag = minlag + 1 : subcount
    cov_ = sum((ts(diag : diag + subseqlen - 1) - mu(diag)) .* (ts(1 : subseqlen) - mu(1)));
    for row = 1 : subcount - diag + 1
        col = diag + row - 1;
        if row > 1
            cov_ = cov_ - dr_bwd(row-1) * dc_bwd(col-1) + dr_fwd(row-1) * dc_fwd(col-1);
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

[motifsIdx] = findMotifs(ts, mu, invnorm, matrixProfile, matrixProfileIdx, subseqlen, minlag);


% Max caps anything that rounds inappropriately. This can hide potential
% issues at times, but it's fairly rare.
matrixProfile = sqrt(max(0, 2 * subseqlen * (1 - matrixProfile), 'includenan'));

% expand initial
if first ~= 1
    isvalidwindow = [false(first - 1, 1); isvalidwindow];
    matrixProfile = [repmat(-1, first - 1, 1); matrixProfile];
    % offset index to compensate for the dropped leading windows
    matrixProfileIdx = [repmat(-1, first - 1, 1); matrixProfileIdx + first - 1];
end

if last ~= initial_subcount
    extendby = initial_subcount - last + 1;
    isvalidwindow = [isvalidwindow; false(extendby, 1)];
    matrixProfile = [matrixProfile; NaN(extendby, 1)];
    matrixProfileIdx = [matrixProfileIdx; NaN(extendby, 1);];
end

% if plot_output
%     gui = mpgui_ed(timeseries, subseqlen);
%     % first pair is plottted alongside data
%     gui.plotProfile(matrixProfile);
%     if isfinite(motifsIdx(1, 1)) && ~isnan(motifsIdx(2, 1))
%         gui.plotData(motifsIdx(1, 1), motifsIdx(2, 1));
%         motiftitle = sprintf('The first motif pair is located at %d and %d.', motifsIdx(1,1), motifsIdx(2,1));
%         gui.plotMotif(gui.motifAx1, motifsIdx(1, 1), motifsIdx(2:end, 1), motiftitle);
%         if isfinite(motifsIdx(1, 2))
%             motiftitle = sprintf('The second motif pair is located at %d and %d.', motifsIdx(1,2), motifsIdx(2,2));
%             gui.plotMotif(gui.motifAx2, motifsIdx(1, 2), motifsIdx(2:end, 2), motiftitle);
%         end
%         if isfinite(motifsIdx(1, 3))
%             motiftitle = sprintf('The third motif pair is located at %d and %d.', motifsIdx(1,3), motifsIdx(2,3));
%             gui.plotMotif(gui.motifAx3, motifsIdx(1, 3), motifsIdx(2:end, 3), motiftitle);
%         end
%         gui.plotDiscords(discordsIdx, sprintf(['The top three discords ', '%d(blue), %d(red), %d(green)'], discordsIdx(1), discordsIdx(2), discordsIdx(3)));
%     else
%         gui.plotData();
%     end
%
%     gui.drawgui;
%
% end

if transposed_  % matches the profile and profile index but not the motif or discord index to the input format
    matrixProfile = transpose(matrixProfile);
    matrixProfileIdx = transpose(matrixProfileIdx);
end

end

function [timeseries, isvalidwindow, first, last] = find_valid_windows_SelfJoin(timeseries, subseqlen)

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



function [mp_a, mp_b, mpi_a, mpi_b] = mpxABBA(a, b, w)
% (mpx_ABBA_v2 titled from Kaveh Kamgar, 2020-12-24 16:56)
% Code and update formulas are by Kaveh Kamgar.
% GUI and top k motif critera are based on but not identical to some code
% by Michael Yeh. Implementation details are specified above.
%
% Additional References
% Yan Zhu, et al, Matrix Profile II: Exploiting a Novel Algorithm and GPUs to break the one Hundred Million Barrier for Time Series Motifs and Join
% Zachary Zimmerman, et al, Scaling Time Series Motif Discovery with GPUs: Breaking the Quintillion Pairwise Comparisons a Day Barrier. (pending review)
% Philippe Pebay, et al, Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments
% Takeshi Ogita, et al, Accurate Sum and Dot Product


% difference equations have 0 as their first entry here to simplify index
% calculations slightly. Alternatively, it's also possible to swap this to the last element
% and reorder the comparison step (or omit on the last step). This is a
% special case when comparing a single time series to itself. The more general
% case with time series A,B can be computed using difference equations for
% each time series.

if isrow(a)
    a = transpose(a);
elseif ~iscolumn(a)
    error('Expected a 1D input for input array a');
end

if isrow(b)
    b = transpose(b);
elseif ~iscolumn(b)
    error('Expected a 1D input for input array b');
end

% matrix profile using cross correlation,
% updated to use the function muinvn, which is a nicer implementation of my
% take on Ogita's work

% save original lengths
subcount_a = length(a) - w + 1;
subcount_b = length(b) - w + 1;

% slightly misleading name. This both identifies valid windows and
% perturbs flat data with a bit of noise, scaled appropriately. If it flips
% any nearest neighbors, they were already unreliable matches.
[a, isvalid_a, first_a, last_a] = find_valid_windows(a, w);
[b, isvalid_b, first_b, last_b] = find_valid_windows(b, w);

mu_a = moving_mean(a, w);
mu_b = moving_mean(b, w);

invn_a = NaN(length(a) - w + 1, 1);
invn_b = NaN(length(b) - w + 1, 1);

% These are truncated so that they start and end on valid subsequences.
% Without that, any adjustment for missing data is susceptible to invalid
% or perturbed initial comparisons
for i = 1 : length(invn_a)
    if isvalid_a(i)
        invn_a(i) = 1 / norm(a(i : i + w - 1) - mu_a(i));
    end
end

for i = 1 : length(invn_b)
    if isvalid_b(i)
        invn_b(i) = 1 / norm(b(i : i + w - 1) - mu_b(i));
    end
end

% This means that for some subsequence, 1/(norm(subseq - mean(subseq))
% is not a finite value, AND this subsequence was not noticed by a mising
% data or constant window check. It is intended as an indicator that you
% should investigate the data itself, because we can't concretely identify
% the issue here. It can be something as simple as intermediate overflow,
% but with double precision inputs, even that is typically indicative
% of a problem with the input.

if any(~isfinite(invn_a) & isvalid_a)
    warning('Times series a has non-normalizable subsequences, which were not detected by the missing data passes');
end

if any(~isfinite(invn_b) & isvalid_b)
    warning('Times series a has non-normalizable subsequences, which were not detected by the missing data passes');
end

invn_a(~isvalid_a) = NaN;
invn_b(~isvalid_b) = NaN;

df = [0; (1/2)*(a(1 + w : end) - a(1 : end - w))];
dg = [0; (a(1 + w : end) - mu_a(2 : end)) + (a(1 : end - w) - mu_a(1 : end - 1))];
dx = [0; (1/2)*(b(1 + w : end) - b(1 : end - w))];
dy = [0; (b(1 + w : end) - mu_b(2 : end)) + (b(1 : end - w) - mu_b(1 : end - 1))];

mp_a = repmat(-1, length(a) - w + 1, 1);
mp_b = repmat(-1, length(b) - w + 1, 1);
mp_a(~isvalid_a) = NaN;
mp_b(~isvalid_b) = NaN;
mpi_a = NaN(length(a) - w + 1, 1);
mpi_b = NaN(length(b) - w + 1, 1);

[mp_a, mp_b, mpi_a, mpi_b] = run_mpx(a, b, mu_a, mu_b, invn_a, invn_b, df, dg, dx, dy, mp_a, mp_b, mpi_a, mpi_b, w);
[mp_b, mp_a, mpi_b, mpi_a] = run_mpx(b, a, mu_b, mu_a, invn_b, invn_a, dx, dy, df, dg, mp_b, mp_a, mpi_b, mpi_a, w);


mp_a = sqrt(2 * w * max(0, 1 - mp_a, 'includenan'));
mp_b = sqrt(2 * w * max(0, 1 - mp_b, 'includenan'));

if first_a ~= 1
    mp_a = [NaN(first_a - 1, 1); mp_a];
    mpi_a = [NaN(first_a - 1, 1); mpi_a];
    % offset indices
    mpi_b = mpi_b + (first_a - 1);
end

if first_b ~= 1
    mp_b = [NaN(first_b - 1, 1); mp_b];
    mpi_b = [NaN(first_b - 1, 1); mpi_b];
    % offset indices
    mpi_a = mpi_a + (first_b - 1);
end

if last_a ~= subcount_a
    mp_a = [mp_a; NaN(subcount_a - length(mp_a), 1)];
    mpi_a = [mpi_a; NaN(subcount_a - length(mp_a), 1)];
end

if last_b ~= subcount_b
    mp_b = [mp_b; NaN(subcount_b - length(mp_b), 1)];
    mpi_b = [mpi_b; NaN(subcount_b - length(mp_b), 1)];
end


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

    c = 1/64;
    % c is just to ensure this doesn't perturb leading digits
    % in case this data provides a partial window for some normalizable
    % section. Feel free to make it smaller or larger.
    % The data variance dependent portion is handled by "scale".
    % 1/64 just happens to be a change in exponent only for any number
    % that doesn't hit the de-normal range. I used it for something close
    % to .01 to avoid perturbing leading digits for the most part.
    timeseries(i : j + subseqlen - 1) = timeseries(i : j + subseqlen - 1) + c * scale * randn(j + subseqlen - i, 1);

    % If any kind of error happens here, this is the point where we give
    % up, as it becomes very difficult to achieve predictable results.
    if any(~isfinite(timeseries(i : j + subseqlen - 1)))
        error(fprintf('unable to fix missing spanning %d to %d\n', i, j));
    end

    i = j + 1;

end

end


function [mpa, mpb, mpia, mpib] = run_mpx(a, b, mua, mub, invna, invnb, df, dg, dx, dy, mpa, mpb, mpia, mpib, w)
amx = length(a) - w + 1;
bmx = length(b) - w + 1;

for ia = 1 : amx
    mx = min(amx - ia + 1, bmx);
    c = sum((a(ia : ia + w - 1) - mua(ia)) .* (b(1 : w) - mub(1)));
    for ib = 1 : mx
        c = c + df(ib + ia - 1) * dy(ib) + dg(ib + ia - 1) * dx(ib);
        c_cmp = c * invna(ib + ia - 1) * invnb(ib);
        if c_cmp > mpa(ib + ia - 1)
            mpa(ib + ia - 1) = c_cmp;
            mpia(ib + ia - 1) = ib;
        end
        if c_cmp > mpb(ib)
            mpb(ib) = c_cmp;
            mpib(ib) = ib + ia - 1;
        end
    end
end
end
