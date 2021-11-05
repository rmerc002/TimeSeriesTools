function [motifIdxs, matrixProfile] = findMotifs(timeSeries, mu, invnorm, matrixProfile, profileIndex, subseqLen, motifCount, neighborCount, exclusionLen, radius)
% This is adapted match the output of some inline code written by Michael Yeh
% to find the top k motifs in a time series. Due to some bug fixes I applied, the two
% may return slightly different results, although they almost always agree on the top 2 motif pairs. 


% +2 allows us to pack the original motif pair with its neighbors.
% Annoyingly if it's extended, matlab will pad with zeros, leading to some
% rather interesting issues later.
motifIdxs = NaN(neighborCount + 2, motifCount);
padLen = 2^nextpow2(length(timeSeries));
crosscov = @(idx) ifft(fft(timeSeries, padLen) .* conj(fft((timeSeries(idx : idx + subseqLen - 1) - mu(idx)) .* invnorm(idx), padLen)), 'symmetric');


for i = 1 : motifCount
    [corr_, motIdx] = max(matrixProfile);
    % -1 means this is maximally negatively correlated and was therefore
    % never updated
    if ~isfinite(corr_) || corr_ == -1 
        break;
    end
    % order subsequence motif pair as [time series index of 1st appearance, time series index of 2nd appearance]
    motifIdxs(1 : 2, i) = [min(motIdx, profileIndex(motIdx)), max(motIdx, profileIndex(motIdx))];
    [corrProfile] = crosscov(motIdx);
    corrProfile = min(1, corrProfile(1 : length(timeSeries) - subseqLen + 1) .* invnorm, 'includenan');
    % This uses correlation instead of normalized Euclidean distance, because it's easier to work with and involves fewer operations.
    
    corrProfile(isnan(matrixProfile)) = NaN;
    if exclusionLen > 0
        for j = 1 : 2
            exclRangeBegin = max(1, motifIdxs(j, i) - exclusionLen + 1); 
            exclRangeEnd = min(length(matrixProfile), motifIdxs(j, i) + exclusionLen - 1);
            corrProfile(exclRangeBegin : exclRangeEnd) = NaN;
        end
    end
    
    for j = 3 : neighborCount + 2
        [neighborCorr, neighbor] = max(corrProfile);
        if  ~isfinite(neighborCorr) || neighborCorr <= 0 || (1 - neighborCorr >= radius * (1 - corr_)) 
            break;
        end
        motifIdxs(j, i) = neighbor;
        if exclusionLen > 0
            exclRangeBegin = max(1, neighbor - exclusionLen + 1); 
            exclRangeEnd = min(length(matrixProfile), neighbor + exclusionLen - 1);
            corrProfile(exclRangeBegin : exclRangeEnd) = NaN;
            corrProfile(exclRangeBegin : exclRangeEnd) = NaN;
        end
    end
    % Matlab allows NaNs to be ignored from min/max calculations, so this
    % is a Matlab specific way to iteratively add excluded ranges of elements. 
    matrixProfile(isnan(corrProfile)) = NaN;
end
end