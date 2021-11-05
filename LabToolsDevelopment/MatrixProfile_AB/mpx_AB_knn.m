function [mp, mpi] = mpx_AB_knn(a, b, w, k)

n = length(a);

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
mp = repmat(-1, k,n);
mpi = NaN(k,n);

for ia = 1 : amx
    fprintf("ia: %d of %d\n",ia, amx);
    mx = min(amx - ia + 1, bmx);
    c = sum((a(ia : ia + w - 1) - mua(ia)) .* (b(1 : w) - mub(1)));
    for ib = 1 : mx
        c = c + df(ib + ia - 1) * dy(ib) + dg(ib + ia - 1) * dx(ib);
        c_cmp = c * invna(ib + ia - 1) * invnb(ib);
        if c_cmp > min(mp(:,ib + ia - 1))
            [mp(:,ib + ia - 1), mpi(:,ib + ia - 1)] = discardMin(c_cmp, ib, mp(:,ib + ia - 1),mpi(:,ib + ia - 1), ceil(w/2));
        end
    end
end

for ib = 1 : bmx
    fprintf("ib: %d of %d\n",ib, bmx);
    mx = min(bmx - ib + 1, amx);
    c = sum((b(ib : ib + w - 1) - mub(ib)) .* (a(1 : w) - mua(1)));
    for ia = 1 : mx
        c = c + df(ia) * dy(ib + ia - 1) + dg(ia) * dx(ib + ia - 1);
        c_cmp = c * invna(ia) * invnb(ib + ia - 1);
        if c_cmp > min(mp(:,ia))
            [mp(:,ia), mpi(:,ia)] = discardMin(c_cmp, ia + ib - 1, mp(:,ia),mpi(:,ia), ceil(w/2));
        end
    end
end
if(any(mp > 1))
    warning('possible precision loss due to rounding');
end
% mp = sqrt(2 * w * (1 - min(1, mp, 'includenan')));
% mp = max(0,sqrt(1 - min(1, mp, 'includenan')));
mp = abs(1-mp);
% mp = max(0,sqrt(1-mp));


end

function [mp,mpi] = discardMin(dist, index, mp, mpi, minlag)
    %insert new distance into existing vector
    if sum(dist > mp) > 0
        %need to make sure the top k are at least minlag apart
        indexProximities = abs(mpi-index);
        [minIndexProximity, discardIndex] = min(indexProximities);
        if ~isempty(discardIndex) && minIndexProximity <= minlag
            if dist > mp(discardIndex)
                %result is not sorted
                mp(discardIndex) = dist;
                mpi(discardIndex) = index;
                [mp, sortedIndices] = sort(mp,"descend");
                mpi = mpi(sortedIndices);
            end
        else %just remove the smallest one
            [tempDistances, sortedIndices] = sort([mp;dist],"descend");
            mp = tempDistances(1:end-1);
            tempProfileIndices = [mpi;index];
            mpi = tempProfileIndices(sortedIndices(1:end-1));
        end
    end
end
