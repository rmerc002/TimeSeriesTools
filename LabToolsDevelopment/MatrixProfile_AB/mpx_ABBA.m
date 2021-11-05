function [mpa, mpb, mpia, mpib] = mpx_ABBA(a, b, w)

if size(a,1) ~= numel(a)
    if size(a,2) == numel(a)
        a = transpose(a);
    else
        error('multidimensional inputs are not supported');
    end
end
if size(b,1) ~= numel(b)
    if size(b,2) == numel(b)
        b = transpose(b);
    else
        error('multidimensional inputs are not supported');
    end
end

% for now assume no nans or consts

% matrix profile using cross correlation,
% updated to use the function muinvn, which is a nicer implementation of my
% take on Ogita's work


nanmapA = movsum(isnan(a), [0 w-1], 'Endpoints', 'discard') ~= 0;
nanmapB = movsum(isnan(b), [0 w-1], 'Endpoints', 'discard') ~= 0;

% This avoids having to constantly break calculations. We just remember
% which subsequences contain > 0 NaNs
a(isnan(a)) = 0;
b(isnan(b)) = 0;


[mua, invna] = muinvn(a,w);
[mub, invnb] = muinvn(b,w);


% this will ensure vectors of constants are effectively uncorrelated with
% everything else. The intention is to maintain consistency with the
% behavior of zscore(,1) when used on a vector of constants, not with the
% behavior of corr(vector_of_constants, vector), which would be NaN
invna(isinf(invna)) = 0;
invnb(isinf(invnb)) = 0;

% differentials have 0 as their first entry. This simplifies index
% calculations slightly and allows us to avoid special "first line"
% handling.

df = [0; (1/2)*(a(1 + w : end) - a(1 : end - w))];
dg = [0; (a(1 + w : end) - mua(2 : end)) + (a(1 : end - w) - mua(1 : end - 1))];
dx = [0; (1/2)*(b(1 + w : end) - b(1 : end - w))];
dy = [0; (b(1 + w : end) - mub(2 : end)) + (b(1 : end - w) - mub(1 : end - 1))];

mpa = repmat(-1, length(a) - w + 1, 1);
mpb = repmat(-1, length(b) - w + 1, 1);
mpia = NaN(length(a) - w + 1, 1);
mpib = NaN(length(b) - w + 1, 1);

% It might be possible to simplify this further in order to avoid
% potential subtle bugs. We want to ensure that we can update smoothly,
% yet anything that touches NaN remains NaN.

if any(nanmapA) || any(nanmapB)
    mpa(nanmapA) = NaN;
    mpb(nanmapB) = NaN;
    [mpa, mpb, mpia, mpib] = solver(a, b, mua, mub, invna, invnb, df, dg, dx, dy, mpa, mpb, mpia, mpib, nanmapA, nanmapB, w);
    [mpb, mpa, mpib, mpia] = solver(b, a, mub, mua, invnb, invna, dx, dy, df, dg, mpb, mpa, mpib, mpia, nanmapB, nanmapA, w);
    mpa = sqrt(max(0, 1 - mpa));
    mpa(nanmapA) = NaN;
    mpb = sqrt(max(0, 1 - mpb));
    mpb(nanmapB) = NaN;
else
    [mpa, mpb, mpia, mpib] = solver_simple(a, b, mua, mub, invna, invnb, df, dg, dx, dy, mpa, mpb, mpia, mpib, w);
    [mpb, mpa, mpib, mpia] = solver_simple(b, a, mub, mua, invnb, invna, dx, dy, df, dg, mpb, mpa, mpib, mpia, w);
    mpa = sqrt(max(0, 1 - mpa));
    mpb = sqrt(max(0, 1 - mpb));
end


end

function [mpa, mpb, mpia, mpib] = solver_simple(a, b, mua, mub, invna, invnb, df, dg, dx, dy, mpa, mpb, mpia, mpib, w)
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


function [mpa, mpb, mpia, mpib] = solver(a, b, mua, mub, invna, invnb, df, dg, dx, dy, mpa, mpb, mpia, mpib, nanmapa, nanmapb, w)
amx = length(a) - w + 1;
bmx = length(b) - w + 1;

for ia = 1 : amx
    mx = min(amx - ia + 1, bmx);
    c = sum((a(ia : ia + w - 1) - mua(ia)) .* (b(1 : w) - mub(1)));
    for ib = 1 : mx
        c = c + df(ib + ia - 1) * dy(ib) + dg(ib + ia - 1) * dx(ib);
        if nanmapa(ib + ia - 1) || nanmapb(ib)
            continue;
        end
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
