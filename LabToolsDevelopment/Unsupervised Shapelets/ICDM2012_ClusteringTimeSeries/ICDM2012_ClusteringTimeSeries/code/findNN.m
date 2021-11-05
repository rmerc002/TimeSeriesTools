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