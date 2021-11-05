function dm = mpx_distanceMatrix(a,minlag,w)
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
dm = repmat(-1,n-w+1,n-w+1);


for diag = minlag+1:diagmax
    c = (sum((a(diag:diag+w-1)-mu(diag)).*(a(1:w)-mu(1))));
    for offset = 1:n-w-diag+2
        c = c + df(offset)*dg(offset+diag-1) + df(offset+diag-1)*dg(offset);
        c_cmp = c*(sig(offset)*sig(offset+diag-1));

        dm(offset,offset+diag-1) = c_cmp;

        dm(offset+diag-1,offset) = c_cmp;
    end
end
% to do ed

dm = min(1,1-dm);
dm = dm.*~eye(size(dm));

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
