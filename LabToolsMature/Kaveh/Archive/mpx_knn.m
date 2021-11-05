function [mp,mpi] = mpx_knn(a,minlag,w,k)
% a is time series
% minlag is the exclusion zone
% w is the window size
% matrix profile using cross correlation, 

% depends on files sum2s, musigtest, dot2s
n = length(a);
[mu, sig] = muinvn(a,w);

a = reshape(a,length(a),1);

% differentials have 0 as their first entry. This simplifies index
% calculations slightly and allows us to avoid special "first line"
% handling.
df = [0; (1/2)*(a(1+w:n)-a(1:n-w))];
dg = [0; (a(1+w:n) - mu(2:n-w+1)) + (a(1:n-w)-mu(1:n-w))];
diagmax = length(a)-w+1;
mp = repmat(-1, k, n);
mpi = NaN(k, n);

for diag = minlag+1:diagmax
    fprintf("%d of %d\n",diag, diagmax);
    c = (sum((a(diag:diag+w-1)-mu(diag)).*(a(1:w)-mu(1))));
    for offset = 1:n-w-diag+2
        c = c + df(offset)*dg(offset+diag-1) + df(offset+diag-1)*dg(offset);
        c_cmp = c*(sig(offset)*sig(offset+diag-1));
        
        if c_cmp > min(mp(:, offset))
            [mp(:, offset), mpi(:, offset)] = discardMin(c_cmp, offset+diag-1, mp(:,offset), mpi(:,offset), minlag);
        end
        
        if c_cmp > min(mp(:, offset+diag-1))
            [mp(:, offset+diag-1), mpi(:, offset+diag-1)] = discardMin(c_cmp, offset, mp(:,offset+diag-1), mpi(:,offset+diag-1), minlag);
        end
    end
end
% to do ed
mp = abs(1-mp);

% mp = sqrt(2*w*(1-mp));
% mp = mp/sqrt(2*w);

% mp = 2*w*(1-mp)/(w*log(w));

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