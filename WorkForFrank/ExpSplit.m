%%
% Returns a binary tree traversal of the integers from 1 to n
function [IDX] = ExpSplit(Range)
n = length(Range);
% Preallocate memory to store the traversed indices and the remaining intervals
IDX       = [1];
% tempPowers = [0];

lowPower = 0;
highPower = log(n);
numSamples = 2;
% powerRange = highPower-lowPower;
while length(IDX) < n
    powers = linspace(lowPower,highPower,numSamples);
    indices = unique(round(exp(powers)));
    uniqueIndices = indices(~ismember(indices,IDX));
    if length(uniqueIndices) == 0
        break;
    end
    IDX = [IDX,uniqueIndices];
    numSamples = numSamples * 2;
%     powerStep = powerRange/(numSamples+1);
%     tempN = length(IDX);
%     for i = 1:tempN
%         newPower = tempPowers(i)+powerStep;
%         newIdx = ceil(exp(newPower));
%         if ~ismember(newIdx,IDX) && newIdx <= Range(end)
%             tempPowers = [tempPowers, newPower];
%             IDX = [IDX, newIdx];
%         end
%     end
%     numSamples = numSamples * 2;
end
IDX(2:end) = IDX(end:-1:2);

end
