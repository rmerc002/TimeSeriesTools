function [cr] = KRASS(a, q, invn)
% The name is a joke.
% precompute fft(ts, 2^nextpow2(length(ts))); and pass that as parameter 1
% pass in (q - mu(q)) .* invn(q) for the second parameter
% 

cv = ifft(fft(a,2^nextpow2(length(a))) .* conj(fft(q, 2^nextpow2(length(a)))), 'symmetric');
[cr] = cv(1 : end - length(q) + 1) .* invn;
cr = sqrt(2*length(q)*(1-cr));
end

