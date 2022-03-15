% constructDemoData;
[simmat] = SimMat(pos, mm); % self-join

%% Remove the diagnomal which will be trivial matches
exclusionLength = ceil(mm/2);
for rr = 1:size(simmat,1)
    startIndex = max(1,rr-exclusionLength+1);
    endIndex = min(size(simmat,1), rr+exclusionLength-1);
    simmat(startIndex:endIndex,rr) = 2*sqrt(mm);
end

distPercentile = 50;
%%% start with a small block size and increase until blocks are solid
blockSize = 10;
simmatBlockMin = matrixPool(simmat, blockSize, "min");

simmatBlockNorm = thresholdMatrix(simmatBlockMin, distPercentile, "min");
figure;
imagesc(simmatBlockNorm);
colormap("gray");



%%% Contrastive Step
MP_AB = mpx_ABBA_v2(pos, neg, mm);

simmatContrast = zeros(size(simmat));
for rr = 1:size(simmat,1)
    for cc = 1:size(simmat,2)
        minAB = min(MP_AB(rr), MP_AB(cc));
        normAB = min(sqrt(2*mm), minAB);
        normAA = min(sqrt(2*mm), simmat(rr,cc));
        simmatContrast(rr,cc) = max(0, normAB - normAA);
    end
end


distPercentile = 10;
%%% start with a small block size and increase until blocks are solid
blockSize = 10;
simmatBlockMin = matrixPool(simmatContrast, blockSize, "max");

simmatBlockNorm = thresholdMatrix(simmatBlockMin, distPercentile,'max');
figure;
imagesc(1-simmatBlockNorm);
colormap("gray");




