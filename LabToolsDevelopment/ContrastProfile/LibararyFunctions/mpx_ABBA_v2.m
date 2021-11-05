function [mp_a, mp_b, mpi_a, mpi_b] = mpx_ABBA_v2(a, b, w)

% Code and update formulas are by Kaveh Kamgar. 
% GUI and top k motif critera are based on but not identical to some code
% by Michael Yeh. Implementation details are specified above.
%
% Additional References
% Yan Zhu, et al, Matrix Profile II: Exploiting a Novel Algorithm and GPUs to break the one Hundred Million Barrier for Time Series Motifs and Join
% Zachary Zimmerman, et al, Scaling Time Series Motif Discovery with GPUs: Breaking the Quintillion Pairwise Comparisons a Day Barrier. (pending review)
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
