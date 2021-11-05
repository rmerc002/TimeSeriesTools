function [plato, top1candidate, top2candidate] = PLATO(positiveTS, negativeTS, subLength, groundTruthIndices)
    %%%Ryan Mercer
    %%%Test Version for Eamonn
    %%%
    %%%Requirements: statistics and machine learning (zscore)
    %%%
    %%% Change Log:
    %%% 2020-12-30, 17:39
    %%% mpx functions now support NaN values
    %%% magnitude conversion taken place outside of mpx
    %%%
    %%% 2020-12-21, 14:45
    %%% Plato is now the top1, not a derived subsequence
    %%% Acceptes input as row and column vector
    %%% 
    %%% 2020-12-16, 14:08
    %%% Matrix profiles use euclidean distances, clipped values above
    %%%    sqrt(2*subLength)
    %%% Distance conversion from matrix profiles to contrast profile
    %%%   divide difference by sqrt(2*subLength)
    %%%   ignore values below 0
    showPlot = true;
    
    
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
    
    

    [MP_AA, MP_AA_Indices] = mpx_plato(positiveTS,subLength,subLength);
    MP_AA = real(MP_AA);
    [MP_AB, MP_AB_Indices] = mpx_AB_plato(positiveTS, negativeTS, subLength);
    MP_AB = real(MP_AB);

    minLength = min([length(MP_AA),length(MP_AB)]); %TODO, this is fishy, I think lengths should be the same
    ContrastProfile = MP_AB(1:minLength) - MP_AA(1:minLength);
    ContrastProfile = ContrastProfile/(sqrt(2*subLength)); %normalize between 0,1
    %discard negative values. These occur when a subsequence in T+ has a
    %closer nearest neighbor in T-. It's not a behavior of interest for
    %this research.
    ContrastProfile = max(0, ContrastProfile); 
    
    [maxContrastValue, maxContrastIndex] = max(ContrastProfile);
    
    fprintf('Top1: %d\n', maxContrastIndex);

    plato = positiveTS(maxContrastIndex:maxContrastIndex + subLength - 1);
    if dataOrientation == 1
       plato = plato'; 
    end

    %%%%%%%%%%%%%
    %%% PLOTS %%%
    %%%%%%%%%%%%%
    if showPlot == true
        redColor = [0.73,0.05,0];
        greenColor = [0,0.73,0.41]; 
        blueColor = [0,0.29,0.73];
        grayColor = [0.65,0.65,0.65];
        lightBlueColor = [0.01, 0.83,0.99];

        tsLength = length(positiveTS);
        maxTSLength = max(length(positiveTS),length(negativeTS));
        %%% for debugging
        fig = figure; 
        set(gcf, 'Position', [0,100,2050,600]);
        
        tiledlayout(6,1);
        
        ax1 = nexttile;
%         subplot(5,1,1);
        plot(negativeTS);
        formattedTitle = sprintf("Negative Time Series, Length=%d",length(negativeTS));
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        xlim([1,maxTSLength]);

        ax2 = nexttile;
%         subplot(5,1,2);
        plot(positiveTS);
        hold on;
        for i = 1:length(groundTruthIndices)
            gti = groundTruthIndices(i);
            %%%TODO: Change the size back to subLength instead of
            %%%2*subLength
            %%% when ground truth is not concatenated
            plot(gti:gti+subLength-1,positiveTS(gti:gti+subLength-1),'Color',greenColor);
        end
        hold off;
        xlim([1,maxTSLength]);
        formattedTitle = sprintf("Positive Time Series, Length=%d",length(positiveTS));
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        xlim([1,maxTSLength]);
        
        ax3 = nexttile;
%         subplot(5,1,3);
        plot(MP_AA,'Color',blueColor);
        hold on;
        plot(MP_AB,'Color',redColor);
        xlim([1,maxTSLength]);
        formattedTitle = "MP_AB(red), MP_AA(blue)";
        title(formattedTitle,'Interpreter','none');
        set(gca,'xtick',[],'ytick',[]);

        ax4 = nexttile;
%         subplot(5,1,4);
        plot(ContrastProfile,'Color',grayColor);
        hold on;
        scatter(maxContrastIndex, maxContrastValue,10,'MarkerFaceColor',redColor,'MarkerEdgeColor',redColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 
        hold off;
        xlim([1,maxTSLength]);
        formattedTitle = sprintf("Contrast Profile");
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        
        ax5 = nexttile;
%         subplot(5,1,5);
        plot(0,0);
        hold on;
        for gti = groundTruthIndices
            plot([gti,gti],[0,1],'Color',greenColor);
        end
        xlim([1,tsLength]);
        formattedTitle = sprintf("GroundTruth");
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        
        ax6 = nexttile;
        distanceProfile = mpx_AB_plato(positiveTS, plato, subLength);
        distanceProfile = real(distanceProfile);
        plot(distanceProfile);
        xlim([1,tsLength]);
        formattedTitle = sprintf("Distance Profile: positiveTS join plato");
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        
        linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x')
        
  
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Plot Effectiveness of Candidate Selection %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        numCandidatesToShow = 1;
        x = zeros(numCandidatesToShow,2*subLength);
        y = zeros(size(x));

        cmappingGlobal = false(1,length(positiveTS));
        for i = 1:length(groundTruthIndices)
            gti = groundTruthIndices(i);
            %%%TODO: Change the size back to subLength instead of
            %%%2*subLength
            %%% when ground truth is not concatenated
            cmappingGlobal(gti:gti+subLength-1) = true;
        end
        

        %%%%%%%%%%%%%%%%%%
        %%% PLATO plot %%%
        %%%%%%%%%%%%%%%%%%
        fig = figure;
        set(gcf, 'Position', [0,100,600,200]);

        plot(plato);

        xlim([0,subLength]);
        formattedTitle = sprintf("PLATO, Length=%d",subLength);
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);

    end
end


function [mp,mpi] = mpx_plato(a,minlag,w)
% a is time series
% minlag is the exclusion zone
% w is the window size
% matrix profile using cross correlation,

% depends on files sum2s, musigtest, dot2s
n = length(a);
[mu, sig] = muinvn(a,w);


% differentials have 0 as their first entry. This simplifies index
% calculations slightly and allows us to avoid special "first line"
% handling.
df = [0; (1/2)*(a(1+w:n)-a(1:n-w))];
dg = [0; (a(1+w:n) - mu(2:n-w+1)) + (a(1:n-w)-mu(1:n-w))];
diagmax = length(a)-w+1;
mp = repmat(-1,n-w+1,1);
mpi = NaN(n-w+1,1);

for diag = minlag+1:diagmax
    c = (sum((a(diag:diag+w-1)-mu(diag)).*(a(1:w)-mu(1))));
    for offset = 1:n-w-diag+2
        c = c + df(offset)*dg(offset+diag-1) + df(offset+diag-1)*dg(offset);
        c_cmp = c*(sig(offset)*sig(offset+diag-1));
        if c_cmp > mp(offset)

            mp(offset) = c_cmp;
            mpi(offset) = offset+diag-1;
        end
        if c_cmp > mp(offset+diag-1)
            mp(offset+diag-1) = c_cmp;
            mpi(offset+diag-1) = offset;
        end
    end
end
% to do ed
mp = sqrt(2*w*(1-mp));
mp = min(sqrt(2*w),mp); %clip the anti-correlated values
mp = real(mp);
end

% Functions here are based on the work in
% Ogita et al, Accurate Sum and Dot Product

function [mu,sig] = muinvn(a,w)
% results here are a moving average and stable inverse centered norm based
% on Accurate Sum and Dot Product, Ogita et al


mu = sum2s(a,w)./w;
sig = zeros(length(a) - w + 1, 1);

for i = 1:length(mu)
    sig(i) = sq2s(a(i : i + w - 1) - mu(i));
end

sig = 1./sqrt(sig);
end


function res = sq2s(a)
h = a .* a;
c = ((2^27) + 1) * a;  % <-- can be replaced with fma where available
a1 = (c - (c - a));
a2 = a - a1;
a3 = a1 .* a2;
r = a2 .* a2 - (((h - a1 .* a1) - a3) - a3);
p = h(1);
s = r(1);
for i = 2 : length(a)
    x = p + h(i);
    z = x - p;
    s = s + (((p - (x - z)) + (h(i) - z)) + r(i));
    p = x;
end
res = p + s;
end


function [x,y] = TwoSquare(a)
x = a .* a;
c = ((2^27) + 1) .* a;
a1 = (c - (c - a));
a2 = a - a1;
a3 = a1 .* a2;
y = a2 .* a2 - (((x - a1 .* a1) - a3) - a3);
end

function [ res ] = sum2s(a,w)
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
end

function [ res ] = sum2s_v2(a,w)
res = zeros(length(a) - w + 1, 1);
accum = a(1);
resid = 0;
for i = 2 : w
    m = a(i);
    p = accum;
    accum = accum + m;
    q = accum - p;
    resid = resid + ((p - (accum - q)) + (m - q));
end
res(1) = accum + resid;
for i = w + 1 : length(a)
    m = a(i - w);
    n = a(i);
    p = accum - m;
    q = p - accum;
    r = resid + ((accum - (p - q)) - (m + q));
    accum = p + n;
    t = accum - p;
    resid = r + ((p - (accum - t)) + (n - t));
    res(i - w + 1) = accum + resid;
end
end


function [mp, mpi] = mpx_AB_plato(a, b, w)

% matrix profile using cross correlation,
% updated to use the function muinvn, which is a nicer implementation of my
% take on Ogita's work

[mua, invna] = muinvn(a,w);
[mub, invnb] = muinvn(b,w);

% differentials have 0 as their first entry. This simplifies index
% calculations slightly and allows us to avoid special "first line"
% handling.

df = [0; (1/2)*(a(1 + w : end) - a(1 : end - w))];
dg = [0; (a(1 + w : end) - mua(2 : end)) + (a(1 : end - w) - mua(1 : end - 1))];
dx = [0; (1/2)*(b(1 + w : end) - b(1 : end - w))];
dy = [0; (b(1 + w : end) - mub(2 : end)) + (b(1 : end - w) - mub(1 : end - 1))];

amx = length(a) - w + 1;
bmx = length(b) - w + 1;
mp = repmat(-1, amx, 1);
mpi = NaN(amx, 1);

for ia = 1 : amx
    mx = min(amx - ia + 1, bmx);
    c = sum((a(ia : ia + w - 1) - mua(ia)) .* (b(1 : w) - mub(1)));
    for ib = 1 : mx
        c = c + df(ib + ia - 1) * dy(ib) + dg(ib + ia - 1) * dx(ib);
        c_cmp = c * invna(ib + ia - 1) * invnb(ib);
        if c_cmp > mp(ib + ia - 1)
            mp(ib + ia - 1) = c_cmp;
            mpi(ib + ia - 1) = ib;
        end
    end
end

for ib = 1 : bmx
    mx = min(bmx - ib + 1, amx);
    c = sum((b(ib : ib + w - 1) - mub(ib)) .* (a(1 : w) - mua(1)));
    for ia = 1 : mx
        c = c + df(ia) * dy(ib + ia - 1) + dg(ia) * dx(ib + ia - 1);
        c_cmp = c * invna(ia) * invnb(ib + ia - 1);
        if c_cmp > mp(ia)
            mp(ia) = c_cmp;
            mpi(ia) = ia + ib - 1;
        end
    end
end
if(any(mp > 1))
    warning('possible precision loss due to rounding');
end
% mp = sqrt(2 * w * (1 - min(1, mp, 'includenan')));
% mp = max(0,sqrt(1 - min(1, mp, 'includenan')));
% mp = max(0,sqrt(1-mp));


%mp = sqrt(1-max(0, min(1, mp, 'includenan'),'includenan'));
mp = sqrt(2*w*(1-mp));
mp = min(sqrt(2*w),mp); %clip the anti-correlated values
mp = real(mp);
end

function [topIndices] = topSubsequenceIndices_V03(ts, subsequenceLength)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

exclusionLength = ceil(subsequenceLength);
N = 2*floor(length(ts)/(subsequenceLength));

A = zeros(2,length(ts));
A(1,:) = 1:length(ts);
A(2,:) = ts';
A = A';
B = flip(sortrows(A,2),1);

exclusionZone = zeros(1,length(ts));

topIndices = nan(1,N);
index = 1;
iterNum = 1; %for debugging
while index <= N && iterNum < size(B,1)
    trialIndex = B(iterNum,1);
    if exclusionZone(trialIndex) == 0
        topIndices(index) = trialIndex;
        index = index + 1;
        leftOffset = exclusionLength;
        rightOffset = exclusionLength;
        exclusionZone(max(1,trialIndex - leftOffset):min(length(ts), trialIndex + rightOffset)) = 1;
    end
    iterNum = iterNum + 1;
end

topIndices = topIndices(~isnan(topIndices));
end

