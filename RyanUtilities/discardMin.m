function [mp,mpi] = discardMin(dist, index, mp, mpi, minlag)
    fprintf("test function\n");
    %insert new distance into existing vector
    if sum(dist < mp) > 0
        %need to make sure the top k are at least minlag apart
        indexProximity = abs(mpi-index);
        [~, discardIndex] = min(indexProximity)
        if ~isempty(discardIndex) 
            if dist < mp(discardIndex)
                fprintf("found index to discard: %d", discardIndex);
                %result is not sorted
                mp(discardIndex) = dist;
                mpi(discardIndex) = index;
                [mp, sortedIndices] = sort(mp,"descend");
                mpi = mpi(sortedIndices);
            end
        else %just remove the smallest one
            [tempDistances, sortedIndices] = sort([mp,dist],"descend");
            mp = tempDistances(1:end-1);
            tempProfileIndices = [mpi,index];
            mpi = tempProfileIndices(sortedIndices(1:end-1));
        end
    end
end
