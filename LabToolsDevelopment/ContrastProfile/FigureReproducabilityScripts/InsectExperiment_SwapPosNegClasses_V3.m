%%%I think for preproducability, it's better to work with hard coded scripts
%%% Started 2020-12-18, 10:04am

%%% Original experiment detailed in ContrastProfile_Insect_ScoreTopK_RyanMercer_2020-12-17
%%% Part changed from original experiment on 2020-12-17
%%% flip roles of pos and neg time series

T = LabC060314{1:100:end,3}; %voltages
startIndexPos = 14300;
endIndexPos = 15430;
pos = T(startIndexPos:endIndexPos); %corresponding to 'np', non-probing class
positiveClass = 'c';

startIndexNeg = 13140;
endIndexNeg = 14105;
neg = T(startIndexNeg:endIndexNeg); %corresponding to 'c', stylet passage through plant cells
subLength = 40;

plato = PLATO(pos, neg, subLength, [1 1]);

labels = LabC060314{1:100:end,4} == positiveClass; %binary labels of positive class
%overwrite the training data
Ttest = T;
rng(1); %required for distance profile consistency
Ttest(startIndexPos:endIndexPos) = rand(endIndexPos-startIndexPos+1,1);
Ttest(startIndexNeg:endIndexNeg) = rand(endIndexNeg - startIndexNeg + 1, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get Motif As Query %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
[mp,mpi] = mpx_plato(pos, subLength, subLength);
 
[minVal, minIndex] = min(mp);
motifA = pos(minIndex:minIndex+subLength -1);
minIndexNN = mpi(minIndex);
motifB = pos(minIndexNN:minIndexNN+subLength -1);
motif = (motifA + motifB)/2;

%%% PLOT MOTIF
subsequences = [motif'; motifA'; motifB'];
 
fig = figure;
set(gcf, 'Position', [0,100,600,200]);
% subplot(10,1,[1,9]);
plot(0,0); hold on;
inset = 0.9;
pi = 0;
for i = 1:3
 
 
    if i == 1
        lineWidth = 5;
        plot([1,subLength],[-0.05,-0.05],'--','Color',[0.8,0.8,0.8]);
    else
        lineWidth = 1;
    end
    
    tempTS = subsequences(i,:);
    tempMin = min(tempTS);
    tempMax = max(tempTS);
    tempRange = max(1e-5, tempMax-tempMin);
 
    plot(-pi+inset*(tempTS-tempMin)/tempRange,'LineWidth',lineWidth);
 
    pi = pi + 1;
end
hold off;
xlim([1,subLength]);
formattedTitle = sprintf("Motif Average and Motif Pair. \nLength=%d",subLength);
title(formattedTitle);
set(gca,'xtick',[],'ytick',[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Work with generalized query code %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
queries = [plato'; motif'];
queryLabels = ["PLATO","Motif"];

for queryIndex = 1:size(queries,1)
    query = queries(queryIndex,:)';
    queryLabel = queryLabels(queryIndex);

    dp = mpx_AB_plato(Ttest,query,subLength);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Select Top K Indices %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K = 1000;
    [bestIndices] = KLowestDistanceIndices_V05(dp, subLength, K, subLength);
    bestIndices6 = bestIndices;
    if K ~= length(bestIndices)%sometimes K is too large
        error("Choice of K = %d was too large. %d bestIndices returned.\n",K,length(bestIndices));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Distance Profile Plot %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    tiledlayout(3,1);

    ax1 = nexttile();
    plot(Ttest);
    title("T Test, with training section set to rand");

    ax2 = nexttile();
    plot(labels);
    title("Labels, Class Positive = 'np'");
    ylim([-0.1, 1.1]);

    ax3 = nexttile();
    plot(dp);
    titleString = sprintf("Distance Profile T Test, %s", queryLabel);
    title(titleString);

    linkaxes([ax1 ax2 ax3],'x');



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                        Top K Bin Breakdown                              
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%
    %%% Print Results %%%
    %%%%%%%%%%%%%%%%%%%%%
    %based on binSize = 10;
    top10 = 100*sum(labels(bestIndices(1:10)))/10;
    top100 = 100*sum(labels(bestIndices(1:100)))/100;
    top1000 = 100*sum(labels(bestIndices(1:1000)))/1000;
    defaultRate = sum(labels)/length(labels);
    defaultClass = defaultRate > 0.5; %means positive class is default
    if defaultRate < 0.5
       defaultRate = 1-defaultRate; 
    end

    resultString = "";
    resultString = resultString + sprintf("%s Results\n", queryLabel);
    resultString = resultString + sprintf("%.1f%%: Top K=10\n",top10);
    resultString = resultString + sprintf("%.1f%%: Top K=100\n",top100);
    resultString = resultString + sprintf("%.1f%%: Top K=1000\n",top1000);
    resultString = resultString + sprintf("%.1f%%: Class %d Default Rate\n",100*defaultRate, defaultClass);
    resultString = resultString + sprintf("%.1f%%: Class %d Rate\n",100*(1-defaultRate), 1-defaultClass);
    resultString = resultString + newline();
    disp(resultString);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot Bin Breakdown %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    plotBinBreakdownOfTopK(K, 10, labels, bestIndices, defaultRate, queryLabel);
    plotBinBreakdownOfTopK(K, 100, labels, bestIndices, defaultRate, queryLabel);
end


function plotBinBreakdownOfTopK(K, binSize, labels, bestIndices, defaultRate, queryLabel)
    if mod(K,10) ~= 0
        error("Code assumes Top K is a multiple of 10, K=%d",K);
    end

    x_index = zeros(ceil(length(bestIndices)/binSize),1);
    y_scorePerGroup = zeros(ceil(length(bestIndices)/binSize),1);
    for i=1:binSize:length(bestIndices)
        groupIndex = ceil(i/binSize);
        x_index(groupIndex) = i;
        labelIndices = bestIndices(i:i+binSize-1);
        y_scorePerGroup(groupIndex) = sum(labels(labelIndices));
    end

    figure; 
    bar(x_index,y_scorePerGroup);
    hold on;
    y = binSize*defaultRate;
    plot([0,K],[binSize*defaultRate, binSize*defaultRate],'--r');
    hold off;
    titleString = sprintf("Bin Breakdown of Top %d, binSize = %d, using %s as query",K,binSize, queryLabel);
    title(titleString);
end
