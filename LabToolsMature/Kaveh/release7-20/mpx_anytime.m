function [matrixProfile, profileIndex, motifIndex, discordIndex] = mpx_anytime(T, exclusion, motif_length)

% Inputs are a time series T, an exclusion window length, and a motif
% length

% exclusion affects the following things
% First, for matrixProfile(i), i - profileIndex(i) >= exclusion
% If position k contains a neighbor or discord, k - exclusion + 1 ..... k + exclusion - 1 
% does not contain a neighbor of any motif or a discord. 

% Outputs are a profile, index, and the locations of each motif and discord


transposed = false;
if ~isvector(T)
    error('malformatted time series');
elseif isrow(T)
    T = transpose(T);
    transposed = true;
end

fprintf('Your time series is of length %d.\n', length(T));

softlgminlen = 13;

% We need to set a downsampling bound based on the log of the motif length
% and the length of the overall time series. We just choose whatever is
% stricter.
upperdsFactor = max(floor(log2(motif_length)) - 3, 0);
upperdsFactor = min(upperdsFactor, max(0, min(softlgminlen, ceil(log2(length(T))) - softlgminlen)));

if upperdsFactor == 0
    disp('This time series or motif length is too short to warrant an anytime algorithm, we are computing it directly.. ');
    t = tic();
    [matrixProfile, profileIndex, mu, invnorm] = mpx(T, exclusion, motif_length);
    elapsed = toc(t);
    [motifIndex, discordIndex] = getMotifsDiscords(T, mu, invnorm, matrixProfile, profileIndex, exclusion, motif_length);
    fprintf('Native complete, time taken:%g seconds\n', elapsed);
else
    stride = 2^upperdsFactor;
    motlenaugm = ceil(motif_length / stride);
    exclusionaugm = ceil(exclusion / stride);
    Taugm = T(1 : stride : end);
    fprintf('\nWe have downsampled the length of your time series to %d or about 1 in %d.\nWe have also shortened the requested motif length to %d and the exclusion length to %d to match.\n\n', length(Taugm), stride, motlenaugm, exclusionaugm);
    t = tic();
    [matrixProfile, profileIndex, mu, invnorm] = mpx(Taugm, exclusionaugm, motlenaugm);
    elapsed = toc(t);
    [motifIndex, discordIndex, updateHandle] = getMotifsDiscords(Taugm, mu, invnorm, matrixProfile, profileIndex, exclusionaugm, motlenaugm);
    fprintf('Length %d Complete, This took time:%g seconds. Press the gui stop button to accept this result.\n', length(Taugm), elapsed);
    
    for i = upperdsFactor - 1 : -1 : 0
        drawnow;
        shouldHalt = updateHandle.shouldHalt();
        if shouldHalt
            break;
        end
        t = tic();
        % mu and invnorm indicate the sliding window mean and sliding
        % window norm(subseq - mean(subseq)). Since these are only used
        % internally, they only match the final matrixProfile and its index on the last iteration
        % if the remaining output if mpx is not preempted.
        if i ~= 0
            stride = 2^i;
            motlenaugm = ceil(motif_length / stride);
            exclusionaugm = ceil(exclusion / stride);
            Taugm = T(1 : stride : end);
            fprintf('\nWe have downsampled the length of your time series to %d or about 1 in %d.\nWe have also shortened the requested motif length to %d and the exclusion length to %d to match.\n\n', length(Taugm), stride, motlenaugm, exclusionaugm);
            [matrixProfileupd, profileIndexupd, mu, invnorm, preempted] = mpx(Taugm, exclusionaugm, motlenaugm, updateHandle);
        else
            fprintf('\nWe are now computing at native length.\n');
            [matrixProfileupd, profileIndexupd, mu, invnorm, preempted] = mpx(T, exclusion, motif_length, updateHandle);
        end
        elapsed = toc(t);
        if ~preempted
            updateHandle.close();
            matrixProfile = matrixProfileupd;
            profileIndex = profileIndexupd;
            if i ~= 0
                [motifIndex, discordIndex, updateHandle] = getMotifsDiscords(Taugm, mu, invnorm, matrixProfile, profileIndex, exclusionaugm, motlenaugm);
                guiTitle = sprintf('1 in %d', stride);
                updateHandle.setTitle(guiTitle);
                fprintf('Currently showing the resampled result for length: %d. This took %g seconds to complete. Press the stop button on the gui to accept this result\n', length(Taugm), elapsed);
            else
                [motifIndex, discordIndex, updateHandle] = getMotifsDiscords(T, mu, invnorm, matrixProfile, profileIndex, exclusion, motif_length);
                fprintf('Currently showing the result for native length: %d. This took %g seconds to complete.\n', length(T), elapsed);
                updateHandle.setTitle('Native');
            end
            drawnow;
        end
    end
    if transposed
        matrixProfile = transpose(matrixProfile);
        profileIndex = transpose(profileIndex);
    end
end
matrixProfile = sqrt(2 * motif_length * (1 - min(1, matrixProfile, 'includenan')));
end

function [motifIndex, discordIndex, guiHandle] = getMotifsDiscords(timeSeries, mu, invnorm, matrixProfile, profileIndex, exclusion, motif_length)
[motifIndex, mpAugmented] = findMotifs(timeSeries, mu, invnorm, matrixProfile, profileIndex, motif_length, 3, 10, exclusion, 2);
[discordIndex] = findDiscords(mpAugmented, 3, exclusion);
matrixProfile = sqrt(2 * motif_length * (1 - min(1, matrixProfile, 'includenan')));
guiHandle = mpgui.launchGui(timeSeries, mu, invnorm, matrixProfile, motifIndex, discordIndex, motif_length);
drawnow;
end
