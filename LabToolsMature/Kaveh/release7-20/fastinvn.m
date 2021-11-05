function [invn] = fastinvn(ts, mu, sublen)
% This is a simple variation on Welford's method.
% This version still results in some minor cancellation and could be
% improved. It isn't prone to anything disastrous.
invn = zeros(length(mu), 1);
invn(1) = sum((ts(1 : sublen) - mu(1)).^2);
for i = 2:length(invn)
    invn(i) = invn(i - 1) + ((ts(i - 1) - mu(i - 1)) + (ts(i + sublen - 1) - mu(i))) * (ts(i + sublen - 1) - ts(i - 1));
end
invn = 1./sqrt(invn);
end