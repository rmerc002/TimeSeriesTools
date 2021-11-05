function [mp, mpi] = mpx_AB_shapeletDiscovery(a, b, w)

% matrix profile using cross correlation,
% updated to use the function muinvn, which is a nicer implementation of my
% take on Ogita's work

[mua, invna] = muinvn(a,w);
[mub, invnb] = muinvn(b,w);

% differentials have 0 as their first entry. This simplifies index
% calculations slightly and allows us to avoid special "first line"
% handling.

df = [0; (1/2)*(a(1 + w : end) - a(1 : end - w))];
dg = [0; (a(1 + w : end) - mua(2 : end)) + (a(1 : end - w) - mua(1 : end - 1))];
dx = [0; (1/2)*(b(1 + w : end) - b(1 : end - w))];
dy = [0; (b(1 + w : end) - mub(2 : end)) + (b(1 : end - w) - mub(1 : end - 1))];

amx = length(a) - w + 1;
bmx = length(b) - w + 1;
mp = repmat(-1, amx, 1);
mpi = NaN(amx, 1);

for ia = 1 : amx
    mx = min(amx - ia + 1, bmx);
    c = sum((a(ia : ia + w - 1) - mua(ia)) .* (b(1 : w) - mub(1)));
    for ib = 1 : mx
        c = c + df(ib + ia - 1) * dy(ib) + dg(ib + ia - 1) * dx(ib);
        c_cmp = c * invna(ib + ia - 1) * invnb(ib);
        if c_cmp > mp(ib + ia - 1)
            mp(ib + ia - 1) = c_cmp;
            mpi(ib + ia - 1) = ib;
        end
    end
end

for ib = 1 : bmx
    mx = min(bmx - ib + 1, amx);
    c = sum((b(ib : ib + w - 1) - mub(ib)) .* (a(1 : w) - mua(1)));
    for ia = 1 : mx
        c = c + df(ia) * dy(ib + ia - 1) + dg(ia) * dx(ib + ia - 1);
        c_cmp = c * invna(ia) * invnb(ib + ia - 1);
        if c_cmp > mp(ia)
            mp(ia) = c_cmp;
            mpi(ia) = ia + ib - 1;
        end
    end
end
if(any(mp > 1))
    warning('possible precision loss due to rounding');
end
% mp = sqrt(2 * w * (1 - min(1, mp, 'includenan')));
% mp = max(0,sqrt(1 - min(1, mp, 'includenan')));
% mp = max(0,sqrt(1-mp));

mp = 1-max(0, min(1, mp, 'includenan'),'includenan');
end

