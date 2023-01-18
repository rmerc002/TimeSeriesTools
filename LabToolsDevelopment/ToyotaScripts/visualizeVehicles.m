% dataConcatAnom = loadVehicleData("Data_anom");
% dataConcatBase = loadVehicleData("Data_base");


%%% 5: Speed of Anomalous Vehicle
%%% 6: Distance traveled
%%% 7: 
% columnNumber = 5;
% upsampleRate = 1;%5;
% 
% posnan = dataConcatAnom(:,columnNumber);
% % posnan = resample(posnan,upsampleRate,1);
% posnan = posnan + randn(length(posnan),1)*1e-5;
% 
% 
% negnan = dataConcatBase(:,columnNumber);
% % negnan = resample(negnan,upsampleRate,1);
% negnan = negnan + randn(length(negnan),1)*1e-5;
% 
% sampleRate = 10;
% mm = 5*sampleRate;
% threshold = 0.1;
% ep = OnlineEmergenceProfile(posnan, negnan, mm, threshold);
% ep.plot();
% ep.plotNoveletDistances();

%%% Distance between following and anomalous
upsampleRate = 1;%5;
X1 = dataConcatAnom(:,3);
X2 = dataConcatAnom(:,8);
Y1 = dataConcatAnom(:,4);
Y2 = dataConcatAnom(:,9);
posnan = sqrt((X1-X2).^2 + (Y1-Y2).^2);
% posnan = resample(posnan,upsampleRate,1);
posnan = posnan + randn(length(posnan),1)*1e-5;

X1 = dataConcatBase(:,3);
X2 = dataConcatBase(:,8);
Y1 = dataConcatBase(:,4);
Y2 = dataConcatBase(:,9);
negnan = sqrt((X1-X2).^2 + (Y1-Y2).^2);
% negnan = resample(negnan,upsampleRate,1);
negnan = negnan + randn(length(negnan),1)*1e-5;

% sampleRate = 10;
% mm = 5*sampleRate;
% threshold = 0.25;
% ep = OnlineEmergenceProfile(posnan, negnan, mm, threshold);
% ep.plot();
% ep.plotNoveletDistances();