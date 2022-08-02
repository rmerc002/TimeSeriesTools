function TimeSeriesChains (A, SubsequenceLength, anchor)
%%Anonymous-Author information blinded for review
%%This is the source code for the ICDM paper "Time Series Chains: A New Primitive
%%for Time Series Data Mining". For more details, please refer to the
%%supporting website: https://sites.google.com/site/timeserieschain/

%%input:
%%A: Time Series
%%SubsequenceLength: Subsequence Length
%%anchor: location of the anchor subsequence. If this input is empty, we
%%will output the unanchored time series chain.

%%output: 
%%Two figures, one showing the chain within the time series, the other
%%enumerating all the subsequence in the chain, and at the same time
%%showing the difference between every two subsequences in the chain.
%%

%% set trivial match exclusion zone
exclusionZone = round(SubsequenceLength/4);

%% check input
if SubsequenceLength > length(A)/2
    error('Error: Time series is too short relative to desired subsequence length');
end
if SubsequenceLength < 4
    error('Error: Subsequence length must be at least 4');
end
if length(A) == size(A, 2)
   A = A'; 
end

%% initialization
MatrixProfileLength = length(A) - SubsequenceLength + 1;
MPLeft = repmat(inf,MatrixProfileLength, 1);
MPindexLeft = zeros(MatrixProfileLength, 1);
MPRight = repmat(inf,MatrixProfileLength, 1);
MPindexRight = zeros(MatrixProfileLength, 1);
[X, n, sumx2, sumx, meanx, sigmax2, sigmax] = ...
    fastfindNNPre(A, SubsequenceLength);

%% compute the matrix profile
pickedIdx = 1:MatrixProfileLength; 
dropval=0;
distanceProfile=zeros(MatrixProfileLength,1);
lastz=zeros(MatrixProfileLength,1);
updatePosLeft=false(MatrixProfileLength,1);
updatePosRight=false(MatrixProfileLength,1);

for i = 1:MatrixProfileLength
    % compute the distance profile
    
    idx = pickedIdx(i);
    subsequence = A(idx:idx+SubsequenceLength-1);
    
    
    if i==1
        [distanceProfile(:,1) lastz dropval lastsumy lastsumy2]= fastfindNN(X, subsequence, n, SubsequenceLength, ...
            sumx2, sumx, meanx, sigmax2, sigmax);
        distanceProfile(:,1) = abs(distanceProfile);
        firstz=lastz;
    else
        lastsumy = lastsumy-dropval+subsequence(end);
        lastsumy2 = lastsumy2-dropval^2+subsequence(end)^2;
        meany=lastsumy/SubsequenceLength;
        sigmay2 = lastsumy2/SubsequenceLength-meany^2;
        sigmay = sqrt(sigmay2);
        lastz(2:n-SubsequenceLength+1)=lastz(1:n-SubsequenceLength)-A(1:n-SubsequenceLength)*dropval+A(SubsequenceLength+1:n)*subsequence(SubsequenceLength);
        %lastz(1)=sum(A(1:SubsequenceLength).*subsequence);
        lastz(1)=firstz(i);
        distanceProfile(:,1) = sqrt(2*(SubsequenceLength-(lastz-SubsequenceLength*meanx*meany)./(sigmax*sigmay)));
        %distanceProfile = sqrt(distanceProfile);
        dropval=subsequence(1);
    end
    
    % apply exclusion zone
    exclusionZoneStart = max(1, idx-exclusionZone);
    exclusionZoneEnd = min(MatrixProfileLength, idx+exclusionZone);
    distanceProfile(exclusionZoneStart:exclusionZoneEnd) = inf;
    
    % figure out and store the neareest neighbor
    updatePosLeft(1:(idx-1)) = false;
    updatePosLeft(idx:MatrixProfileLength) = distanceProfile(idx:MatrixProfileLength) < MPLeft(idx:MatrixProfileLength);
    MPLeft(updatePosLeft) = distanceProfile(updatePosLeft);
    MPindexLeft(updatePosLeft) = idx;
    
    updatePosRight((idx+1):MatrixProfileLength) = false;
    updatePosRight(1:idx) = distanceProfile(1:idx) < MPRight(1:idx);
    MPRight(updatePosRight) = distanceProfile(updatePosRight);
    MPindexRight(updatePosRight) = idx;
    
    
    if mod(i,1000)==0
        i
    end
    
end

ChainPos = false(MatrixProfileLength,1);
ChainLength = zeros(MatrixProfileLength,1);

for i = 1:(MatrixProfileLength-1)
    if (~ChainPos(i))
        cur=i;
        count=1;
        while MPindexRight(cur)>0 && MPindexLeft(MPindexRight(cur))==cur
            ChainPos(cur)=true;
            cur=MPindexRight(cur);
            count=count+1;
        end
        ChainLength(i)=count;
    end
end

%%Note that ChainLength and ChainPos show all possible chains
%%within the time series. The following code only outputs the longest chain
%%possible, but actually with these two vectors, we can show any chain
%%within the time series. 

ChainStart=0;
if (~exist('anchor'))
    ChainStart=find(ChainLength==max(ChainLength),1);
    %ChainStart=ChainStart(1);
else
    ChainStart=anchor;
end

count=ChainLength(ChainStart);

figure1=figure();
plot(A);
hold on;

figure2=figure();

curmax=0;
curmin=0;
curmaxdiff=0;
curmindiff=0;

cur=ChainStart;
next=ChainStart;
lastpattern=[];
for i=1:count
    next=MPindexRight(cur);
    curpattern=zscore(A(cur:(cur+SubsequenceLength-1)),1);
    if max(curpattern)>curmax
        curmax=max(curpattern);
    end
    if min(curpattern)<curmin
        curmin=min(curpattern);
    end
    if i>=2
        diff=curpattern-lastpattern;
        val1=max(diff);
        val2=min(diff);
        if val1>curmaxdiff
            curmaxdiff=val1;
        end
        if val2<curmindiff
            curmindiff=val2;
        end
    end
    lastpattern=curpattern;
    cur=next;
end

cur=ChainStart;
i=1;
lastpattern=[];

while i<=count
    figure(figure2);
    subplot(ceil(count/3),3,i);
    curpattern=zscore(A(cur:(cur+SubsequenceLength-1)),1);
    plot(cur:(cur+SubsequenceLength-1),curpattern);
    xlim([cur,cur+SubsequenceLength-1]);
    ylim([curmin,curmax]);
    if i>1
        hold on;
        plot(cur:(cur+SubsequenceLength-1),-lastpattern+curpattern,'r');
    end
    figure(figure1);
    plot(cur:(cur+SubsequenceLength-1),A(cur:(cur+SubsequenceLength-1)),'r');
    cur=MPindexRight(cur);
    lastpattern=curpattern;
    i=i+1;
end

%figure;
%cur=ChainStart;
    
%for i=1:count-1
%    next=MPindexRight(cur);
%    subplot(ceil(count/3),3,i);
%    plot(cur:(cur+SubsequenceLength-1),-zscore(A(cur:(cur+SubsequenceLength-1)),1)+zscore(A(next:(next+SubsequenceLength-1)),1));
%    xlim([cur,cur+SubsequenceLength-1]);
%    ylim([curmindiff,curmaxdiff]);
%    if i==1
%        title('evolution trend');
%    end
%    cur=next;
%end

% m is winSize
function [X, n, sumx2, sumx, meanx, sigmax2, sigmax] = fastfindNNPre(x, m)
n = length(x);
x(n+1:2*n) = 0;
X = fft(x);
cum_sumx = cumsum(x);
cum_sumx2 =  cumsum(x.^2);
sumx2 = cum_sumx2(m:n)-[0;cum_sumx2(1:n-m)];
sumx = cum_sumx(m:n)-[0;cum_sumx(1:n-m)];
meanx = sumx./m;
sigmax2 = (sumx2./m)-(meanx.^2);
sigmax = sqrt(sigmax2);

% m is winSieze
function [dist lastz dropval sumy sumy2] = fastfindNN(X, y, n, m, sumx2, sumx, meanx, sigmax2, sigmax)
%x is the data, y is the query
%y = (y-mean(y))./std(y,1);                      %Normalize the query
dropval=y(1);
y = y(end:-1:1);                                %Reverse the query
y(m+1:2*n) = 0;

%The main trick of getting dot products in O(n log n) time
Y = fft(y);
Z = X.*Y;
z = ifft(Z);

%compute y stats -- O(n)
sumy = sum(y);
sumy2 = sum(y.^2);
meany=sumy/m;
sigmay2 = sumy2/m-meany^2;
sigmay = sqrt(sigmay2);


dist = 2*(m-(z(m:n)-m*meanx*meany)./(sigmax*sigmay));
dist = sqrt(dist);
lastz=real(z(m:n));