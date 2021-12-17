function [platos, firstIndices, platoContrast, CP] = LeftContrastProfile_Batch(positiveTS, negativeTS, m, contrastThreshold, forcePlot,txt)
    %%% Input:
    %%%   positiveTS: A time series containing at least two instances of a desired behavior
    %%%   negativeTS: A time series containing zero instances of a desired behavior
    %%%   m: The approximate window size of the desired behavior
    %%%   forcePlot: (optional) Display the following two plots
    %%%
    %%% Output:
    %%%   plato: The subsequence that most distinguishes 
    %%%          positiveTS from negativeTS
    %%%   platoIndices: The starting index of each of the K Platos
    %%%   platoPrimaryContrast: Contrast value of each plato in the K=1 CP
    %%%   platoNaryContrast: Contrast value of each plato after appending
    %%%     the previous Plato to negativeTS. May be helpful in identifying
    %%%     diminishing returns and retundant behaviors.
    
    dataOrientation = 0; %0 for row, 1 for column
    %Change to row vector for internal use
    if size(positiveTS,1) == 1 
        %save the data orientation for matching output format to input
        %give orientation priority to positiveTS
        dataOrientation = 1; 
        positiveTS = positiveTS';
    end
    if size(negativeTS,1) == 1 
        %give orientation priority to positiveTS, do not save negative
        %orientation if different than positive
        negativeTS = negativeTS';
    end
    
    
    %%% Matrix profile self-join using positiveTS
    MP_AA = nan(length(positiveTS),1);
    [~, tempMP_AA, ~, MP_AA_Indices] = mpxLeftRight(positiveTS, ceil(m/2), m);
    tempMP_AA = real(tempMP_AA);
    MP_AA(1:length(tempMP_AA)) = clipMatrixProfileAmplitude(tempMP_AA, m);

    %%% Matrix profile AB-join between positiveTS and negativeTS
    historicMP_AB = nan(length(positiveTS),1);
    tempMP_AB = mpxABBA(positiveTS, negativeTS, m);
    historicMP_AB(1:length(tempMP_AB)) = clipMatrixProfileAmplitude(tempMP_AB, m);
    
    MP_AB = historicMP_AB;
    
    %%% Contrast Profile
    historicCP = MP_AB - MP_AA;
    %%% Normalize values to the range [0,1]
    historicCP = normalizeContrastProfileAmplitude(historicCP, m);

    CP = historicCP; %%%CP will be updated (lowered) by comparing to learned platos

    positiveTSNonNan = positiveTS;
    positiveTSNonNan(isnan(positiveTS)) = nanmean(positiveTS);

    negativeTSNonNan = negativeTS;
    negativeTSNonNan(isnan(negativeTS)) = nanmean(negativeTS);

    contextLength = floor(m/2);
    platosWithContext = [];
    platoIndices = [];
    platoNNIndices = [];
    platoContrast = [];

    singleResetFlag = true;

    %%%debugging vars
    oldNewsAnnotation = true(length(positiveTS),1);
    historicOldNewsAnnotation = oldNewsAnnotation;

    numSubsequences = length(positiveTS)-m+1;
    
    for timeIndex = 1:numSubsequences
        
        subsequence = positiveTS(timeIndex:timeIndex + m - 1);
        if sum(isnan(subsequence)) > 0
            continue;
        end

        DPPlato = inf;
        platoDist = inf;
        numPlatos = size(platosWithContext,1);
        for platoIndex = 1:numPlatos
            platoWithContext = platosWithContext(platoIndex,:)';
%             platoDist = getAlignmentDist(platoWithContext, subsequence);
            platoDist = real(MASS_V2(platoWithContext, subsequence));
            platoDist = min(platoDist);
            platoDist = clipMatrixProfileAmplitude(platoDist, m);
            if singleResetFlag == true || platoIndex < numPlatos || timeIndex > platoIndices(platoIndex) + m
                DPPlato = min(DPPlato, platoDist);
            end
        end
        
        ABDist = min(MP_AB(timeIndex), DPPlato);

        MP_AB(timeIndex) = ABDist;


        CPi = ABDist - MP_AA(timeIndex);
        CPi = normalizeContrastProfileAmplitude(CPi, m);

        if ~isempty(platoIndex) > 0
            historicMP_AB(timeIndex) = MP_AB(timeIndex);
            historicCP(timeIndex) = CP(timeIndex);
            historicOldNewsAnnotation(timeIndex) = oldNewsAnnotation(timeIndex);
        end

        %%% If a past subsequence exceeding the contrast threshold, but
        %%% is close to a subsequence with greater contrast, replace it
        if singleResetFlag == true && ~isempty(platosWithContext) && timeIndex < platoIndices(end) + m && CPi > platoContrast(end)
%             fprintf("Replacing plato index:%d with index:%d\n", platoIndices(end), timeIndex);
            startIndex = platoIndices(end);
            endIndex = timeIndex;
            CP(startIndex:endIndex) = historicCP(startIndex:endIndex);
            MP_AB(startIndex:endIndex) = historicMP_AB(startIndex:endIndex);
            oldNewsAnnotation(startIndex:endIndex) = historicOldNewsAnnotation(startIndex:endIndex);

            platosWithContext(end,:) = [];
            platoIndices(end) = [];
            platoNNIndices(end) = [];
            platoContrast(end) = [];
            
        else
            DPPlato = min(DPPlato, platoDist);
            ABDist = min(MP_AB(timeIndex), DPPlato);
            MP_AB(timeIndex) = ABDist;
            CPi = ABDist - MP_AA(timeIndex);
            CPi = normalizeContrastProfileAmplitude(CPi, m);
        end
        
        %%% Normalize values to the range [0,1]
        CP(timeIndex) = CPi;
        
        
        if CPi > contrastThreshold
            %%%TODO: need to handle nans in the AB dist calc when I adjust the bound of timeIndex
%             startIndex = max(1, timeIndex - m);
%             endIndex = min(numSubsequences, timeIndex + 2*m + 1);
            if ~isempty(platoIndices) && timeIndex < platoIndices(end) + m
                singleResetFlag = false;
            else
                singleResetFlag = true;
            end
            startIndex = timeIndex - m;
            endIndex = timeIndex + 2*m + 1;
            %%%Todo: fix, this is temporary so I don't have to handle nans
            %%%in platos with context
            if startIndex < 1 || endIndex > numSubsequences
                continue;
            end
            subsequenceWithContext = positiveTS(startIndex:endIndex);
            platosWithContext(end+1,:) = subsequenceWithContext;
            platoIndices(end+1) = timeIndex;
            platoNNIndices(end+1) = MP_AA_Indices(timeIndex);
            platoContrast(end+1) = CPi;                
            fprintf("Found new behavior at index %d with contrast value: %.2f\n", platoIndices(end), platoContrast(end));
            
            
            
     
        elseif historicCP(timeIndex) > contrastThreshold %%% for debugging
            oldNewsAnnotation(timeIndex:timeIndex+m-1) = false;
        end
        
        
        
    end 

    firstIndices = platoNNIndices;
%     firstIndices = platoIndices;
    numPlatos = size(platosWithContext,1);
    if numPlatos > 0
        platos = nan(numPlatos,m);
        for platoIndex = 1:numPlatos
            startIndex = platoNNIndices(platoIndex);
            endIndex = startIndex+m-1;
            platos(platoIndex,:) = positiveTS(startIndex:endIndex);
        end

        [firstIndices, sortedIndices] = sort(firstIndices);
        platos = platos(sortedIndices,:);
        platoContrast = platoContrast(sortedIndices);
    else
        platos = [];
    end

    MP_AA_Indices(end+1:end+m-1) = length(positiveTS)-m+2:length(positiveTS);
    MP_AA_Indices(1:m) = 1:m;
    xaxis = 1:length(positiveTS);
    nanIndices = isnan(MP_AA_Indices);
    MP_AA_Indices(nanIndices) = xaxis(nanIndices);
    CPOld = CP;
    CP = zeros(size(CP));
    for ii = 1:length(positiveTS)
        CP(MP_AA_Indices(ii)) = max(CP(MP_AA_Indices(ii)), CPOld(ii));
    end
%     CP(MP_AA_Indices) = CP; %%%Assign CP values to the left nearest neighbor index
%     MP_AA(MP_AA_Indices) = MP_AA;
%     MP_AB(MP_AA_Indices) = MP_AB;
%     CP(1:m) = 0;
%     CP(nanIndices) = 0;
    %%%%%%%%%%%%%
    %%% PLOTS %%%
    %%%%%%%%%%%%%
    if forcePlot == true
        
        redColor = [0.73,0.05,0];
        greenColor = [0,0.73,0.41]; 
        blueColor = [0,0.29,0.73];
        grayColor = [0.75,0.75,0.75];
        lightGrayColor = [0.9, 0.9, 0.9];
        lightBlueColor = [0.01, 0.83,0.99];
        platoColor = [129/255, 51/255, 144/255];
        platoTwinColor = [115/255, 170/255, 43/255];

        tsLength = length(positiveTS);
        maxTSLength = max(length(positiveTS),length(negativeTS));

        fig = figure('Name','Contrast Profile: Panel','NumberTitle','off'); 
        set(gcf, 'Position', [0,100,2000,600]);
        
        tiledlayout(4,1);
        
        ax1 = nexttile;
        plot(negativeTS,'Color',redColor);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}T^{(-)}\\color{black}: Negative Time Series",redColor(1), redColor(2), redColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,length(negativeTS)],'ytick',[], 'TickDir','out');
        xlim([1,maxTSLength]);
        box off;
        
        ax2 = nexttile;
        hold on;
        tempPos = positiveTS;
        tempPos(oldNewsAnnotation) = nan;
        plot(tempPos);

        tempPos = positiveTS;
        tempPos(~oldNewsAnnotation) = nan;
        plot(tempPos,'Color',lightGrayColor,'LineWidth',1);

        hold off;
        xlim([1,maxTSLength]);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}T^{(+)}\\color{black}: Positive Time Series",blueColor(1),blueColor(2),blueColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,length(positiveTS)],'ytick',[], 'TickDir','out');
        xlim([1,maxTSLength]);
        box off;
        
        ax3 = nexttile;
        plot(MP_AA,'Color',blueColor);
        hold on;
        plot(MP_AB,'Color',redColor);
        xlim([1,maxTSLength]);
        ylim([0,sqrt(2*m)]);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}MP^{(+ -)} AB-join, \\color[rgb]{%f,%f,%f}MP^{(+ +)} Self-Join",redColor(1), redColor(2), redColor(3),blueColor(1),blueColor(2),blueColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,length(positiveTS)],'ytick',[0,sqrt(2*m)], 'TickDir','out');
        box off;
        
        ax4 = nexttile;
        plot(CP,'Color',grayColor);
        hold on;
        scatter(platoNNIndices, ones(1,length(platoIndices)),10,'MarkerFaceColor',platoColor,'MarkerEdgeColor',platoColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 
        scatter(platoIndices, 1,10,'MarkerFaceColor',platoTwinColor,'MarkerEdgeColor',platoTwinColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0);

        hold off;
        xlim([1,maxTSLength]);
        ylim([0,1]);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}Contrast Profile, \\color[rgb]{%f,%f,%f}First, \\color[rgb]{%f,%f,%f}Repeat",grayColor(1), grayColor(2), grayColor(3), platoColor(1), platoColor(2), platoColor(3), platoTwinColor(1), platoTwinColor(2), platoTwinColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,length(positiveTS)],'ytick',[0,1], 'TickDir','out');
        box off;
        
        
        
        linkaxes([ax1 ax2 ax3 ax4],'x')
        
        %%%Generating unique colors for classes when there can be many classes
        %%% I will assume no more than 1000 classes
        numColors = 1000;
        colors = lines(numColors);
        
        figure;
        hold on;
        numPlatos = size(platos,1);
        for ki = 1:numPlatos
            tempTS = platos(ki,:);
            tempMin = min(tempTS);
            tempMax = max(tempTS);
            tempRange = tempMax-tempMin;

            plot(-ki + 0.9*(tempTS - tempMin)/tempRange);
        end
        hold off;
        formattedTitle = sprintf("Top-%d Platos", numPlatos);
        title(formattedTitle);
        set(gca,'xtick',[1,m],'ytick',[], 'TickDir','out');
        xlim([0,m]);
        box off;
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

