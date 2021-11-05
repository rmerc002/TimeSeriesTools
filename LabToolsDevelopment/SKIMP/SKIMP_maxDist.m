function MAGICMatrixProfile = SKIMP_maxDist(Dataset, varargin)

%   License to use and modify this code is granted freely without warranty to all, as long as the original author is
%   referenced and attributed as such. The original author maintains the right to be solely associated with this work.
%
%   Programmed and Copyright by Frank Madrid: fmadr002[at]ucr[dot]edu
%   Date: 06/07/2019

%% INPUT VALIDATION
p = inputParser;
paramName     = 'Dataset';
validationFcn = @(x) validateattributes(x, {'numeric'}, {'vector'});
addRequired(p, paramName, validationFcn);

paramName     = 'Range';
defaultVal    = 3:1:length(Dataset) / 4;
validationFcn = @(x) validateattributes(x, {'numeric'}, {'integer', 'positive', 'increasing'});
addOptional(p, paramName, defaultVal, validationFcn);

paramName     = 'SaveImage';
defaultVal    = false;
validationFcn = @(x) validateattributes(x, {'logical'}, {'scalar'});
addParameter(p, paramName, defaultVal, validationFcn);

paramName     = 'ShowImage';
defaultVal    = false;
validationFcn = @(x) validateattributes(x, {'logical'}, {'scalar'});
addParameter(p, paramName, defaultVal, validationFcn);

p.parse(Dataset, varargin{:});
INPUTS = p.Results;

assert(max(INPUTS.Range) < length(Dataset) / 2, '[SKIMP] Error: The maximum subsequence length explored must be less than half the size of the dataset.');
%% INITIALIZE

% Determine the order in which we will explore the subsequence lengths used to run the matrix profile
SplitIDX = BinarySplit(length(INPUTS.Range));



% Stores the intermediate values of each Matrix Profile where MAGICMatrixProfile(i,:) = MatrixProfile with subsequence length SplitIDX(i)
MAGICMatrixProfile = nan(length(INPUTS.Range),length(INPUTS.Dataset));

%RYAN
prevIndices = [1,length(INPUTS.Range)];
M = mpx(INPUTS.Dataset,  floor(INPUTS.Range(1) / 4), INPUTS.Range(1));
MAGICMatrixProfile(INPUTS.Range(1),1:length(M)) = abs(M);
  
M = mpx(INPUTS.Dataset,  floor(INPUTS.Range(end) / 4), INPUTS.Range(end));
MAGICMatrixProfile(INPUTS.Range(end),1:length(M)) = abs(M);

% Stores row indices for image
ImageIDX = ones(length(INPUTS.Range),1);

% Initialize and displays the progress bar which tracks SKIMP's status, may be 'Cancelled' to interrupt the processing of SKIMP
ProgressBar = waitbar(0, 'Calculating the MAGIC Profile', 'Name', 'SKIMP', 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(ProgressBar, 'canceling', 0);

if(INPUTS.SaveImage)
  FILENAME = char(randi([97,122],1,64)) % Purposely not suppressed
end

if(INPUTS.ShowImage || INPUTS.SaveImage)
  ViewSkimp = figure('Name', 'MAGIC Matrix Profile', 'Position', [24 300 1800 140]);
  TempAx = axes('Position', [0 0 1 1], 'Parent', ViewSkimp); hold on;
  xlim([1 length(INPUTS.Dataset) - max(INPUTS.Range)]);
  xticks([]);
  yticks([]);
end

if(INPUTS.ShowImage)
  ViewSKIMP.Visible = true;
end
%% BEGIN PROCESSING

% Runs matrix profile on each subsequence length in SplitIDX, stores the results in MAGICMatrixProfile, and generates a frame for the animation
for i = 3:length(INPUTS.Range)
  
  waitbar(i / length(INPUTS.Range), ProgressBar, 'Calculating the MAGIC Profile');
  
  % Premptively terminate if the user cancels the ProgressBar
  if getappdata(ProgressBar,'canceling')
    break
  end
  
  % Get the current subsequence length to run the Matrix Profile on
  metaIndex = getNextIndex(MAGICMatrixProfile, INPUTS.Range, prevIndices);
  if metaIndex == 0
      return;
  end
  disp(metaIndex);
  prevIndices = sort([prevIndices,metaIndex]);
  SubsequenceLength = INPUTS.Range(metaIndex);
  
  
  % Run Matrix Profile
  M = mpx(INPUTS.Dataset,  floor(SubsequenceLength / 4), SubsequenceLength);
  MAGICMatrixProfile(INPUTS.Range(metaIndex),1:length(M)) = abs(M); 
  
  % Propogate changes to ImageIDX
  % TO-DO: Using a binary search to find last index, update all values in range. IDK if this will result in much of a times save
  j = metaIndex;
  upperBoundIndex = prevIndices(find(prevIndices==metaIndex)+1);
  while(j <= length(INPUTS.Range) && j ~= upperBoundIndex)
    ImageIDX(j) = metaIndex;
    j = j + 1;
  end
  % Convert MAGICMatrixProfile to a gray scale image
  [A, map] = rgb2ind(ind2rgb(ceil(size(gray(256),1)*(MAGICMatrixProfile(ImageIDX, :))), gray(256)),256);

  % Draw a red line across the current explored index
  map(end+1,:) = [1 0 0]; % Create a new mapped color
  A(metaIndex,:) = size(map,1)-1;
  
  if(INPUTS.ShowImage || INPUTS.SaveImage)
    set(0, 'Currentfigure', ViewSkimp);
    cla(TempAx);
    imagesc(TempAx, A);
    colormap(map);
    delete(findall(gcf,'type','annotation'))
    annotation('textbox', [0.01, 0.9, 0.1, 0.1], 'String', ['Completion Percentage: ' num2str(i/length(SplitIDX)*100,'%.2f') '%']);
  end
  
  if(INPUTS.SaveImage)
    F = getframe(TempAx);
    Image = frame2im(F);
    [A,map] = rgb2ind(Image,256);
    if i == 1
      imwrite(A, map, [FILENAME '.gif'], 'gif', 'LoopCount', inf, 'DelayTime', 0.001);
    else
      imwrite(A, map, [FILENAME '.gif'], 'gif', 'WriteMode','append','DelayTime',0.001);
    end
  end
    
end
delete(ProgressBar);
end % function SKIMP

function [IDX] = BinarySplit(n)
% BinarySplit(Lowerbound, Upperbound)
% Preallocate memory to store the traversed indices and the remaining intervals
IDX       = [];
Intervals = {};

IDX(1) = 1; % We always begin by explore the first integer
Intervals(1) = {[2,n]}; % After exploring the first integer, we begin splitting the interval 2:n

while(~isempty(Intervals))
    lb = Intervals{1}(1);
    ub = Intervals{1}(2);
    mid = floor((lb + ub) / 2);
    Intervals(1) = [];
    
    IDX(end+1) = mid;
    
    if(lb == ub)
        continue;
    else
        [L,R] = split(lb,ub,mid);
        if(~isempty(L))
            Intervals{end+1} = L;
        end
        if(~isempty(R))
            Intervals{end+1} = R;
        end
    end
end


function [L,R] = split(lb,ub, m)
    if(lb == m)
        L = [];
        R = [m+1, ub];
    elseif ub == m
        L = [lb m-1];
        R = [];
    else
        L = [lb m-1];
        R = [m+1 ub];
    end
end
            
            
end

function [mp,mpi] = mpx(a,minlag,w)

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

function nextIndex = getNextIndex(PMP, range, prevIndices)
    prevIndices = sort(prevIndices);
    maxDist = 0;
    nextIndex = 0;
%     nextIndex = find((prevIndices(2:end)-prevIndices(1:end-1))-1>0);
%     if length(nextIndex) == 0
%         return
%     end
%     nextIndex = nextIndex(1)+1;
    for metaIdx = 2:length(prevIndices)
        idx1 = prevIndices(metaIdx-1);
        idx2 = prevIndices(metaIdx);
        if idx1 + 1 < idx2 %can't search between integers
            thisDist = PMP(range(idx2),:)-PMP(range(idx1),:);
            thisDist = thisDist(~isnan(thisDist));
            thisDist = norm(thisDist)/sqrt(2*length(thisDist));
            if thisDist > maxDist
                maxDist = thisDist;
                nextIndex = ceil((idx2+idx1)/2);
            end
        end
    end
end
