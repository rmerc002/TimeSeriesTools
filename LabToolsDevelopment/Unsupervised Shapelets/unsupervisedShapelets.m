function [feature, location] = unsupervisedShapelets(data)    
    feature = [];
    location = {};
    
    ts = 1;

    data2=data;
    ql=50;

    DIS=[];
    FTR=1;
    iter=1;
    
    while length(feature) == 0
        [r,c]=size(data2);
        cls2 = [];
        cls2(r,1)=0;
        seps=[];
        cnt=1;
        startLength = ceil(0.1*size(data,2));
        endLength = size(data,2);
        numLengths = 50;
        sampleSizes = getSubLenSeries(startLength, endLength, numLengths);
        for ql=sampleSizes   %60,   %sony=10:35
            for i=2:c-ql
               ftr = data2(ts,i:i+ql-1);
               dis = computeMatrix(ftr,data2);
               dis = dis/sqrt(ql);

               seps(cnt,1)=0;   %gap
               seps(cnt,2)=0;   %corr
               seps(cnt,3)=0;   %q=length(ind1)/length(ind2)
               seps(cnt,4)=ql;  %ql= qlen
               seps(cnt,5)=i;   %index
               for corr=0.95:-0.01:0.65        
                    ind1 = find(dis(:,1)<sqrt(2*(1-corr)));
                    ind2 = find(dis(:,1)>sqrt(2*(1-corr)));
                    if(length(ind1)==0 || length(ind2)==0)
                        continue;
                    end
                    m1 = mean(dis(ind1,1)); 
                    s1 = std(dis(ind1,1),1); 
                    m2 = mean(dis(ind2,1));
                    s2 = std(dis(ind2,1),1);

                    q = length(ind1)/length(ind2);

                    if( q > 0.2 && q < 5 )    
                        curr = (m2-s2-(m1+s1));
                    else
                        curr = 0;
                    end
                    if curr>seps(cnt,1)
                        seps(cnt,1) = curr;
                        seps(cnt,2) = corr;
                        seps(cnt,3) = q;
                    end
               end 
               cnt = cnt+1;
            end
            disp(ql);
        end
        [a,b]=sortrows(seps,-1);
        indx = 0;
        for i=1:size(a,1)
            q=a(i,3);
            if(q>0.2 && q<5)  %0.42
                indx = i;
                break;
            end
        end
        if(indx==0)
          break; %indx=1;
        end
        ql=a(indx,4);
        ftr = data2(ts,a(indx,5):a(indx,5)+ql-1);
        dis = computeMatrix(ftr,data2);
        dis = dis/sqrt(ql);
        corr = a(indx,2);
        ind1 = find(dis(:,1)<sqrt(2*(1-corr)));
        ind2 = find(dis(:,1)>sqrt(2*(1-corr)));
        m1 = mean(dis(ind1,1)); 
        s1 = std(dis(ind1,1),1); 
        m2 = mean(dis(ind2,1));
        s2 = std(dis(ind2,1),1);

    %     dis = computeMatrix(ftr,data2); 
        dist=computeMatrix(ftr,data);
        dist=dist/sqrt(ql);
        DIS = [DIS dist];






    %         I = find(dis(:,1)<(m1+0.5*s1));
    %         cls2(I)=1;
    %         X = find(cls2==0);
    %         cls1 = cls1(X);
    %         data2=data2(X,:);
    %         FTR = FTR+1;

        feature = ftr;
        location = {ts,a(indx,5),a(1,2)};



    %         [x,y]=sort(dis,'descend');
    %         for i=1:r
    %             if(cls2(y(i))==0)
    %                 ts=find(X==y(i));
    %                 break;
    %             end
    %         end

    end
end

function dis = computeMatrix(feats, data)
    n = size(data,1);
    F = size(feats,1);
    dis(1:n,1:F) = 0;
    for j = 1:n
        for i = 1:F
            [loc(j,i) dis(j,i)] = findNN(data(j,2:end),feats(i,:));
        end
    end
end

function [loc bsf] = findNN(x,y)
    %x is the data, y is the query
    format long;
    xt=x;
    yt=y;
    n = length(x);
    y = (y-mean(y))./std(y,1); %Normalize the query
    
    m = length(y);
    x(n+1:2*n) = 0;
    y = y(end:-1:1);                                %Reverse the query
    y(m+1:2*n) = 0;
    
    %The main trick of getting dot products in O(n log n) time
    X = fft(x);
    Y = fft(y);
    Z = X.*Y;
    z = ifft(Z);

    %compute y stats -- O(n)
    sumy = sum(y);
    sumy2 = sum(y.^2);
    
    %compute x stats -- O(n)
    
    
    cum_sumx = cumsum(x);
    cum_sumx2 =  cumsum(x.^2);
    sumx2 = cum_sumx2(m+1:n)-cum_sumx2(1:n-m);
    sumx = cum_sumx(m+1:n)-cum_sumx(1:n-m);
    meanx = sumx./m;
    sigmax2 = (sumx2./m)-(meanx.^2);
    sigmax = sqrt(sigmax2);

    %computing the distances -- O(n) time
    dist = (sumx2 - 2*sumx.*meanx + m*(meanx.^2))./sigmax2 - 2*(z(m+1:n) - sumy.*meanx)./sigmax + sumy2;
    dist = sqrt(dist);
    
    %find the minimum
    [bsf index] = min(dist);
    loc = index+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   HELPER FUNCTIONS   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subLenSeries = getSubLenSeries(startLength, endLength, numLengths)
%%% Purpose: Reduce space.assume subsequences are pretty similar from
%%% from one subLen to the next. By trial and error, adding the square root
%%% of the current subLen seems to produce a good distribution.

powerMin = log10(startLength);
powerMax = log10(endLength);
powerStep = (powerMax-powerMin)/numLengths;
powers = powerMin:powerStep:powerMax;
subLenSeries = unique(ceil(power(10,powers)));

end
