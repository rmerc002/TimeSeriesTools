function [score1, score2] = BehaviorSimilarityMeasure(data)
% uiopen('C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\ContrastProfile\Experiments\2021-06-07_PhoneExperiments\Walking_WornNew_MailboxToLightPost\Raw Data.csv',1)
% data = RawData{1:3:end,5};
subLength = 50;
topK = 1 + 1;

figure; plot(data);

endRegime1Index = floor(length(data)/2);

% [mp, mpi] = mpx_AB_knn(data, data, subLength, topK);
[mp, mpi] = mpx_v3(data, subLength);

mpi_regime = mpi >= endRegime1Index;
mpi_regime = mpi_regime*endRegime1Index;
mpi_regime = mpi_regime + 1;

visualizeMMPAB(data, mp(2:end,:), mpi_regime(2:end,:), [], [], [], [subLength], 1:topK-1, "KNN", "KNN");

mpi_class = mpi >= endRegime1Index;

meanWornShoes01 = mean(mpi_class(2:2, 1:endRegime1Index),'all');
meanNewShoes01 = mean(mpi_class(2:2, endRegime1Index:end),'all'); 

% meanWornShoes05 = mean(~mpi_class(2:6, 1:endRegime1Index),'all'); 
% meanNewShoes05 = mean(mpi_class(2:6, endRegime1Index:end),'all'); 
% 
% meanWornShoes10 = mean(~mpi_class(2:11, 1:endRegime1Index),'all'); 
% meanNewShoes10 = mean(mpi_class(2:11, endRegime1Index:end),'all'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Weighted Scores   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
mp_size = size(mp);
mp_flattened = reshape(mp,[numel(mp),1]);
[~,mp_sorted_indices] = sort(mp_flattened);
mp_linspace = linspace(1,0,numel(mp));
mp_uniform_values = nan(numel(mp),1);
mp_uniform_values(mp_sorted_indices) = mp_linspace;
mp_unflattened = reshape(mp_uniform_values, mp_size);

mpi_class = mpi >= endRegime1Index;
mpi_class = mpi_class*2-1;
mpi_class_weighted = mpi_class.*mp_unflattened;

%Top-1 NN
tempValuesNeg = mpi_class_weighted(2:2, 1:endRegime1Index);
tempValuesNeg(tempValuesNeg > 0) = 0;
tempValuesNeg(tempValuesNeg < 0) = -1*tempValuesNeg(tempValuesNeg < 0);
meanWornShoesWeighted01 = mean(tempValuesNeg,'all'); 

tempValuesPos = mpi_class_weighted(2:2, endRegime1Index:end);
tempValuesPos(tempValuesPos < 0) = 0;
meanNewShoesWeighted01 = mean(tempValuesPos,'all'); 

% %Top-5 NN
% tempValuesNeg = mpi_class_weighted(2:6, 1:endRegime1Index);
% tempValuesNeg(tempValuesNeg > 0) = 0;
% tempValuesNeg(tempValuesNeg < 0) = -1*tempValuesNeg(tempValuesNeg < 0);
% meanWornShoesWeighted05 = mean(tempValuesNeg,'all'); 
% 
% tempValuesPos = mpi_class_weighted(2:6, endRegime1Index:end);
% tempValuesPos(tempValuesPos < 0) = 0;
% meanNewShoesWeighted05 = mean(tempValuesPos,'all'); 
% 
% %Top-10 NN
% tempValuesNeg = mpi_class_weighted(2:11, 1:endRegime1Index);
% tempValuesNeg(tempValuesNeg > 0) = 0;
% tempValuesNeg(tempValuesNeg < 0) = -1*tempValuesNeg(tempValuesNeg < 0);
% meanWornShoesWeighted10 = mean(tempValuesNeg,'all'); 
% 
% tempValuesPos = mpi_class_weighted(2:11, endRegime1Index:end);
% tempValuesPos(tempValuesPos < 0) = 0;
% meanNewShoesWeighted10 = mean(tempValuesPos,'all'); 

%Return Values
score1 = meanWornShoesWeighted01;
score2 = meanNewShoesWeighted01;