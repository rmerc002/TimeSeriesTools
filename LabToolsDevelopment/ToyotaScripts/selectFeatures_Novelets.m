columnNumber = 5;

upsampleRate = 1;%5;
posnan = dataConcatAnom(:,columnNumber);
% posnan = resample(posnan,upsampleRate,1);
posnan = posnan + randn(length(posnan),1)*1e-5;

negnan = dataConcatBase(:,columnNumber);
% negnan = resample(negnan,upsampleRate,1);
negnan = negnan + randn(length(negnan),1)*1e-5;

pos = posnan;
pos(isnan(posnan)) = mean(posnan,'omitnan');
% pos = pos + randn(length(pos),1)*1e-5;

neg = negnan;
neg(isnan(negnan)) = mean(negnan,'omitnan');
% neg = neg + randn(length(neg),1)*1e-5;

mm = 10*upsampleRate;
threshold = 0.25;
ep = OnlineEmergenceProfile(posnan, negnan, mm, threshold);
ep.plot();
ep.plotNoveletDistances();

%%% Distance to following vehicle
upsampleRate = 1;%5;
xterm = dataConcatAnom(:,3) - dataConcatAnom(:,8);
yterm = dataConcatAnom(:,4) - dataConcatAnom(:,9);
posnan = sqrt(xterm.^2 + yterm.^2);
posnan = resample(posnan,upsampleRate,1);
posnan = posnan + randn(length(posnan),1)*1e-5;

xterm = dataConcatBase(:,3) - dataConcatBase(:,8);
yterm = dataConcatBase(:,4) - dataConcatBase(:,9);
negnan = sqrt(xterm.^2 + yterm.^2);
negnan = resample(negnan, upsampleRate,1);
negnan = negnan + randn(length(negnan), 1)*1e-5;

pos = posnan;
pos(isnan(posnan)) = mean(posnan,'omitnan');
% pos = pos + randn(length(pos),1)*1e-5;

neg = negnan;
neg(isnan(negnan)) = mean(negnan,'omitnan');
% neg = neg + randn(length(neg),1)*1e-5;

mm = 10*upsampleRate;
threshold = 0.25;
ep = OnlineEmergenceProfile(posnan, negnan, mm, threshold);
ep.plot();