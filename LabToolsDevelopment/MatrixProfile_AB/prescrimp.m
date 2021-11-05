function [mp, mpi] = prescrimp(a, b, w, stride, iters)
%disp(iters);
[mua] = sum2s(a, w) ./ w;
[mub] = sum2s(b, w) ./ w;
[invna] = fastinvn(a, mua, w);
[invnb] = fastinvn(b, mub, w);

malen = length(a) - w + 1;
mblen = length(b) - w + 1;
I = (1 : stride : mblen)';
iters = min(iters, length(I));
I = I(randperm(length(I), iters));
%disp(iters);
fa = fft(a);
mp = repmat(-1, malen, 1);
mpi = zeros(malen, 1);

for i = 1 : iters
    fq = conj(fft(b(I(i) : I(i) + w - 1), length(a))); 
    cv = ifft(fa .* fq, 'symmetric');
    cr = (cv(1 : malen) - w .* mua .* mub(I(i))) .* invna .* invnb(I(i));
    f = find(mp < cr);
    mpi(f) = I(i);
    mp(f) = cr(f);
    [~, mxa] = max(cr);
    mxb = I(i);
    c = cv(mxa);   

    mnj = min(stride, min(mxb - 1, mxa - 1));
    mxj = min(stride, min(malen - mxa, mblen - mxb));
    for j = 1 : mnj 
        c = c - a(mxa - j + w) * b(mxb - j + w) + a(mxa - j) * b(mxb - j); 
        cr = (c - w * mua(mxa - j) * mub(mxb - j)) * invna(mxa - j) * invnb(mxb - j);
        if mp(mxa - j) < cr
            mp(mxa - j) = cr;
            mpi(mxa - j) = mxb - j;
        end
    end
    c = cv(mxa);
    for j = 1 : mxj
        c = c - a(mxa + j - 1) * b(mxb + j - 1) + a(mxa + j + w - 1) * b(mxb + j + w - 1); 
        cr = (c - w * mua(mxa + j) * mub(mxb + j)) * invna(mxa + j) * invnb(mxb + j);
        if mp(mxa + j) < cr
            mp(mxa + j) = cr;
            mpi(mxa + j) = mxb + j;
        end
    end
end
mp = sqrt(max(0, 2 * w * (1 - mp))); % <-- note this is the correct way to avoid complex numbers. abs() is incorrect. 
end

% xsum
function [ res ] = sum2s(p,w)
   % based on kahan-babuska and ACCURATE SUM AND DOT PRODUCT by Ogita et al
   
   mlen = length(p)-w+1;
   res = zeros(mlen,1);
   pi_k = p(1);
   sig_k = 0;
   
   for i = 2:w
       [pi_k,q] = TwoSum(pi_k,p(i));
       sig_k = sig_k + q;
   end
   
   res(1) = pi_k + sig_k;
   
   for i = w+1:length(p)
       pi_k_rem = -1*p(i-w);
       [pi_k,q] = TwoSum(pi_k,pi_k_rem);
       sig_k = sig_k + q;
       [pi_k,q] = TwoSum(pi_k,p(i));
       sig_k = sig_k + q;
       res(i-w+1) = pi_k + sig_k;
   end
   
end

%xadd
function [x,y] = TwoSum(a,b)
   x = a+b; 
   z = x-a;
   y = ((a-(x-z))+(b-z));

end


function [invn] = fastinvn(ts, mu, sublen)
% This is a simple variation on Welford's method. 
% This version still results in some minor cancellation and could be
% improved. It isn't prone to anything disastrous. 
invn = zeros(length(mu),1);
invn(1) = sum((ts(1:sublen)-mu(1)).^2);
for i = 2:length(invn)
    invn(i) = invn(i-1) + ((ts(i-1) - mu(i-1)) + (ts(i+sublen-1) - mu(i))) * (ts(i+sublen-1)-ts(i-1)); 
end
invn = 1./sqrt(invn);
end

