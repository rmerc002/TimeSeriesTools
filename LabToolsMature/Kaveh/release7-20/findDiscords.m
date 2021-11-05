
function [discordIdx] = findDiscords(matrixProfile, discordCount, exclusionLen)
% This function assumes a correlation based implementation of matrixProfile
% and that any excluded regions have been set to NaN, nothing else.
discordIdx = zeros(discordCount, 1);
maxCount = min(2 * discordCount * exclusionLen, numel(~isnan(matrixProfile)));
[corr_, idx] = mink(matrixProfile, maxCount);
for i = 1 : discordCount
    f = find(~isnan(corr_), 1);
    if corr_(f) == -1
        discordIdx(i : end) = NaN;
        break;
    end
    firstValid = idx(f);
    % since we generally initialize to -1, we have to skip here
    if isempty(firstValid) 
        discordIdx(i : end) = NaN;
        break;
    end
    discordIdx(i) = firstValid;
    % this assumes the chosen index is part of the exclusion length, so we
    % have up to exclusionLen - 1 excluded indices on each side. If exclusionLen is
    % greater than or equal to subsequence/motif length, then discords
    % cannot partially overlap.
    corr_(abs(idx - discordIdx(i)) < exclusionLen) = NaN;
end
end