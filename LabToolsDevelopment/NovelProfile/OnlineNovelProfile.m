%%% Input:
%%%   obj.positiveTS: A time series containing at least two instances of a desired behavior
%%%   obj.negativeTS: A time series containing zero instances of a desired behavior
%%%   m: The approximate window size of the desired behavior
%%%   forcePlot: (optional) Display the following two plots
%%%
%%% Output:
%%%   novlet: The subsequence that most distinguishes 
%%%          obj.positiveTS from obj.negativeTS
%%%   obj.novletIndices: The starting index of each of the K Novlets
%%%   novletPrimaryContrast: Contrast value of each novlet in the K=1 obj.NP
%%%   novletNaryContrast: Contrast value of each novlet after appending
%%%     the previous Novlet to obj.negativeTS. May be helpful in identifying
%%%     diminishing returns and retundant behaviors.
classdef OnlineNovelProfile
    properties
        positiveTS
        negativeTS
        mm % Subsequence length
        exclusionLength
        noveltyThreshold
        MP_AA
        MP_AA_Indices
        MP_AB
        CP % Contrast Profile
        NP % Novlet Profile

        minBufferLength
        oldNewsAnnotation
        bufferStartIndex

        novlets
        novletIndices
        novletNNIndices
        novletScores
    end
    methods
        function obj = OnlineNovelProfile(positiveTS, mm, negativeTS, noveltyThreshold)
            if nargin < 1
               positiveTS = [];
            end
            if nargin < 2
                mm = 10;
            end
            if nargin < 3
                negativeTS = [];
            end
            if nargin < 4
                noveltyThreshold = 0.1;
            end
            
            obj.positiveTS = reshape(positiveTS, length(positiveTS), 1);
            obj.mm = mm;
            obj.exclusionLength = obj.mm;
            obj.negativeTS = reshape(negativeTS, length(negativeTS), 1);
            obj.noveltyThreshold = noveltyThreshold;
            
            obj.minBufferLength = obj.mm*3;
            obj.bufferStartIndex = 1;
            obj.oldNewsAnnotation = []; %%% boolean
            obj.novlets = []; %%% dimension: (numNovlets, mm)
            obj.novletIndices = [];
            obj.novletNNIndices = [];
            obj.novletScores = [];

            %%% Matrix profile self-join using obj.positiveTS
            obj.MP_AA = nan(length(obj.positiveTS),1);
            obj.MP_AA_Indices = nan(length(obj.positiveTS),1);
            [~, tempMP_AA, ~, tempMP_AA_Indices] = mpxLeftRight(obj.positiveTS, ceil(obj.mm/2), obj.mm);
            tempMP_AA = real(tempMP_AA);
            obj.MP_AA(1:length(tempMP_AA)) = clipMatrixProfileAmplitude(tempMP_AA, obj.mm);
            obj.MP_AA_Indices(1:length(tempMP_AA_Indices)) = tempMP_AA_Indices;
        
            %%% Matrix profile AB-join between obj.positiveTS and obj.negativeTS
            obj.MP_AB = nan(length(obj.positiveTS),1);
            tempobj.MP_AB = mpxABBA(obj.positiveTS, obj.negativeTS, obj.mm);
            obj.MP_AB(1:length(tempobj.MP_AB)) = clipMatrixProfileAmplitude(tempobj.MP_AB, obj.mm);
            
          

            obj.bufferStartIndex = 1;

            
            obj = discoverNovlets(obj, positiveTS, obj.MP_AB, obj.MP_AA, obj.MP_AA_Indices);
        end

        function obj = discoverNovlets(obj, t, MP_AB, MP_AA, MP_AA_Indices)
            t = reshape(t, length(t), 1);

            
            bufferLength = length(t);
            initialBufferStartIndex = obj.bufferStartIndex;
            while bufferLength >= obj.minBufferLength
                bufferIndex = obj.bufferStartIndex - initialBufferStartIndex + 1;
                startIndexActionWindow = bufferIndex;
                endIndexActionWindow = min(length(obj.positiveTS), startIndexActionWindow + obj.minBufferLength - 1);
                tSub = t(startIndexActionWindow:endIndexActionWindow);
                MP_ABSub = MP_AB(startIndexActionWindow:endIndexActionWindow);
                MP_AASub = MP_AA(startIndexActionWindow:endIndexActionWindow);
                MP_AA_IndicesSub = MP_AA_Indices(startIndexActionWindow:endIndexActionWindow);
                %%% Update the MP_AB here
                for jj = 1:size(obj.novlets,1)
                    novlet = obj.novlets(jj, :)';
                    newMP_AB = real(MASS_V2(tSub, novlet)); %%%TODO: verify if this overlap is correct
                    newMP_AB = clipMatrixProfileAmplitude(newMP_AB, obj.mm); 
                    MP_ABSub(1:length(newMP_AB)) = min(MP_ABSub(1:length(newMP_AB)), newMP_AB);
                end
                
                CPSub = MP_ABSub - MP_AASub;
                CPSub = normalizeContrastProfileAmplitude(CPSub, obj.mm);
                [peakIndices, peakValues] = exclusionZonePeaks(CPSub, obj.mm, obj.exclusionLength);
                peakIndex = peakIndices(1);

                peakMP_AB = MP_ABSub(peakIndex);
                peakMP_AA = MP_AASub(peakIndex);
                peakMP_AA_Index = MP_AA_IndicesSub(peakIndex);

                
                peakContrast = peakMP_AB - peakMP_AA;
                peakContrast = normalizeContrastProfileAmplitude(peakContrast, obj.mm);

                if peakContrast > obj.noveltyThreshold
                    startIndex = peakIndex;
                    endIndex = peakIndex + obj.mm - 1;
                    novlet = tSub(startIndex:endIndex);
                    obj.novlets(end+1,:) = novlet;
                    obj.novletIndices(end+1) = obj.bufferStartIndex + peakIndex - 1;
                    obj.novletNNIndices(end+1) = obj.bufferStartIndex + peakMP_AA_Index - 1;
                    obj.novletScores(end+1) = peakContrast;
                end

                %%%TODO: add NP profile caculated from CP and MP_AA_Indices
                activeEndIndex = obj.bufferStartIndex + length(MP_AASub)-1;
                obj.MP_AA(obj.bufferStartIndex:activeEndIndex) = MP_AASub;
                obj.MP_AB(obj.bufferStartIndex:activeEndIndex) = MP_ABSub;
                obj.CP(obj.bufferStartIndex:activeEndIndex) = CPSub;

                obj.bufferStartIndex = obj.bufferStartIndex + peakIndex + obj.exclusionLength;
                bufferLength = length(t) - (obj.bufferStartIndex - initialBufferStartIndex + 1) + 1;
            end
            

        end

        function obj = update(obj, t)
            if isempty(t)
                return;
            else
                t = reshape(t, length(t), 1);
                obj.positiveTS = [obj.positiveTS; t];
            end
        
            numPossibleSubsequences = length(obj.positiveTS) - obj.mm + 1;
            %%% Assume length(obj.MP_AB) includes m-1 points after last valid value
            currentBufferLength = length(obj.positiveTS) - length(obj.MP_AB);
            if currentBufferLength < obj.minBufferLength
                numSamplesNeededForUpdate = obj.minBufferLength - currentBufferLength;
                fprintf("New samples stored, but %d more samples required to act on buffer", numSamplesNeededForUpdate)
                return;
            end
        
            obj.bufferStartIndex = max(1, obj.bufferStartIndex - 2*obj.exclusionLength);%%%TODO: Temp to debug, but maybe necessary
            unprocessedPositiveTS = obj.positiveTS(obj.bufferStartIndex:end);
            unprocessedLength = length(unprocessedPositiveTS);

            numSubsequences = unprocessedLength - obj.mm + 1;

            newMP_AB = nan(unprocessedLength,1);
            tempMP_AB = mpx_ABBA_v2(unprocessedPositiveTS, obj.negativeTS, obj.mm);
            newMP_AB(1:length(tempMP_AB)) = tempMP_AB;
            %%% Start Novlet Comparisons
            numNovlets = size(obj.novlets,1);
            novletDP = sqrt(2*obj.mm)*ones(numNovlets,unprocessedLength);
            
            if numNovlets > 0
                for novletIndex = 1:numNovlets
                    novlet = reshape(obj.novlets(novletIndex,:), obj.mm, 1);
                    tempNovletDist = real(MASS_V2(unprocessedPositiveTS, novlet));
                    tempNovletDist = clipMatrixProfileAmplitude(tempNovletDist, obj.mm);
                    novletDP(novletIndex,1:numSubsequences) = tempNovletDist;
                end
                [novletMinDist, novletIndices] = min(novletDP,[],1);
                novletMinDist = reshape(novletMinDist, unprocessedLength, 1);
    
                newMP_AB = min(newMP_AB, novletMinDist);
            end


            %%%Extending the left only matrix profile can be done by
            %%%combinine an ab join between unprocessed posTS with rest of
            %%%posTS , then the left only self-join of the unprocessed
            %%%posTS
            
            newMP_AA = sqrt(2*obj.mm)*ones(unprocessedLength, 1);
            newMP_AA_Indices = ones(unprocessedLength, 1);

%             [~, tempLeftAA, ~, tempLeftAAIndices] = mpxLeftRight(obj.positiveTS, ceil(obj.mm/2), obj.mm);
%             tempLeftAA = tempLeftAA(obj.bufferStartIndex:end);
%             newMP_AA(1:length(tempLeftAA)) = tempLeftAA;
            
            [tempAB, ~, tempABIndices] = mpx_ABBA_v2(unprocessedPositiveTS, obj.positiveTS(1:obj.bufferStartIndex-1), obj.mm);
            [~, tempLeftAA, ~, tempLeftAAIndices] = mpxLeftRight(unprocessedPositiveTS, ceil(obj.mm/2), obj.mm);
            tempLeftAAIndices = tempLeftAAIndices + obj.bufferStartIndex -1;
            
            for ii = 1:numSubsequences
                if isnan(tempAB(ii))
                    newMP_AA(ii) = tempLeftAA(ii);
                    newMP_AA_Indices(ii) = tempLeftAAIndices(ii);
                elseif isnan(tempLeftAA)
                    newMP_AA(ii) = tempAB(ii);
                    newMP_AA_Indices(ii) = tempABIndices(ii);
                elseif tempAB(ii) < tempLeftAA(ii) 
                    newMP_AA(ii) = tempAB(ii);
                    newMP_AA_Indices(ii) = tempABIndices(ii);
                else %%% if AB is nan or 
                    newMP_AA(ii) =  tempLeftAA(ii);
                    newMP_AA_Indices(ii) = tempLeftAAIndices(ii);
                end
            end

            obj.bufferStartIndex = obj.bufferStartIndex + 2*obj.exclusionLength;
            obj = discoverNovlets(obj, unprocessedPositiveTS(2*obj.exclusionLength+1:end), newMP_AB(2*obj.exclusionLength+1:end), newMP_AA(2*obj.exclusionLength+1:end), newMP_AA_Indices(2*obj.exclusionLength+1:end));
            %%%end of update
        end
        function plot(obj)
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
                novletColor = [129/255, 51/255, 144/255];
                novletTwinColor = [115/255, 170/255, 43/255];
        
                tsLength = length(obj.positiveTS);
                maxTSLength = max(length(obj.positiveTS),length(obj.negativeTS));
        
                fig = figure('Name','Contrast Profile: Panel','NumberTitle','off'); 
                set(gcf, 'Position', [0,100,2000,600]);
                
                tiledlayout(4,1);
                
                ax1 = nexttile;
                plot(obj.negativeTS,'Color',redColor);
                formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}T^{(-)}\\color{black}: Negative Time Series",redColor(1), redColor(2), redColor(3));
                title(formattedTitle);
                set(gca,'xtick',[1,length(obj.negativeTS)],'ytick',[], 'TickDir','out');
                xlim([1,maxTSLength]);
                box off;
                
                ax2 = nexttile;
                hold on;
                tempPos = obj.positiveTS;
                tempPos(oldNewsAnnotation) = nan;
                plot(tempPos);
        
                tempPos = obj.positiveTS;
                tempPos(~oldNewsAnnotation) = nan;
                plot(tempPos,'Color',lightGrayColor,'LineWidth',1);
        
                hold off;
                xlim([1,maxTSLength]);
                formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}T^{(+)}\\color{black}: Positive Time Series",blueColor(1),blueColor(2),blueColor(3));
                title(formattedTitle);
                set(gca,'xtick',[1,length(obj.positiveTS)],'ytick',[], 'TickDir','out');
                xlim([1,maxTSLength]);
                box off;
                
                ax3 = nexttile;
                plot(obj.MP_AA,'Color',blueColor);
                hold on;
                plot(obj.MP_AB,'Color',redColor);
                xlim([1,maxTSLength]);
                ylim([0,sqrt(2*m)]);
                formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}MP^{(+ -)} AB-join, \\color[rgb]{%f,%f,%f}MP^{(+ +)} Self-Join",redColor(1), redColor(2), redColor(3),blueColor(1),blueColor(2),blueColor(3));
                title(formattedTitle);
                set(gca,'xtick',[1,length(obj.positiveTS)],'ytick',[0,sqrt(2*m)], 'TickDir','out');
                box off;
                
                ax4 = nexttile;
                plot(obj.NP,'Color',grayColor);
                hold on;
                scatter(obj.novletNNIndices, ones(1,length(obj.novletIndices)),10,'MarkerFaceColor',novletColor,'MarkerEdgeColor',novletColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 
                scatter(obj.novletIndices, 1,10,'MarkerFaceColor',novletTwinColor,'MarkerEdgeColor',novletTwinColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0);
        
                hold off;
                xlim([1,maxTSLength]);
                ylim([0,1]);
                formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}Contrast Profile, \\color[rgb]{%f,%f,%f}First, \\color[rgb]{%f,%f,%f}Repeat",grayColor(1), grayColor(2), grayColor(3), novletColor(1), novletColor(2), novletColor(3), novletTwinColor(1), novletTwinColor(2), novletTwinColor(3));
                title(formattedTitle);
                set(gca,'xtick',[1,length(obj.positiveTS)],'ytick',[0,1], 'TickDir','out');
                box off;
                
                
                
                linkaxes([ax1 ax2 ax3 ax4],'x')
                
                %%%Generating unique colors for classes when there can be many classes
                %%% I will assume no more than 1000 classes
                numColors = 1000;
                colors = lines(numColors);
                
                figure;
                hold on;
                numNovlets = size(novlets,1);
                for ki = 1:numNovlets
                    tempTS = novlets(ki,:);
                    tempMin = min(tempTS);
                    tempMax = max(tempTS);
                    tempRange = tempMax-tempMin;
        
                    plot(-ki + 0.9*(tempTS - tempMin)/tempRange);
                end
                hold off;
                formattedTitle = sprintf("Top-%d Novlets", numNovlets);
                title(formattedTitle);
                set(gca,'xtick',[1,m],'ytick',[], 'TickDir','out');
                xlim([0,m]);
                box off;
            end
            %%% end of plot
        end
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
% Zachary Ziobj.mmerman, et al, Scaling Time Series Motif Discovery with GPUs: Breaking the Quintillion Pairwise Comparisons a Day Barrier. (pending review)
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