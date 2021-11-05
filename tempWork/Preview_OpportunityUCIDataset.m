
%plots = zeros(21,floor(plotLength));

data = S1Drill{:,:}';
% data = S1ADL1';
nscale = 1;
plotLength = floor(size(data,2)*nscale);
startIndex = 1000;
endIndex = min(size(data,2)/3, startIndex + plotLength - 1);

cup_accelX_gyroX = [135,138];
fridge_reed = [196, 197, 198];
lowerdrawer_reed = [202, 203, 206];
% chair_accel = [211, 212, 213];
dish_reed = [195, 205, 207];
RLA = [64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76];
back_magnetic = [44, 45, 46];
LRA_accel = [90, 91, 92];
Rshoe_accel = [119, 120, 121];%, 122, 123, 124];
locomotion_label = [244];
highlevel_label = [245];
rightarm_label = [248,249];
leftarm_label = [246,247];
botharms_label = [250];
time = [1];

columns = [];

columns = [columns, cup_accelX_gyroX];
columns = [columns, locomotion_label];
% columns = [columns, highlevel_label];
% columns = [columns, leftarm_label];
columns = [columns, rightarm_label];
columns = [columns, botharms_label];


columns = [columns, dish_reed];
columns = [columns, back_magnetic];

columns = [columns, lowerdrawer_reed];
columns = [columns, RLA];

columns = [columns, fridge_reed];
columns = [columns, LRA_accel];
columns = [columns, Rshoe_accel];


columns = [columns, time];


numPlots = length(columns);

% figure;

plotData = nan(numPlots,plotLength);
for i = 1:numPlots
    dataRow = data(columns(i),startIndex:endIndex);
    dataRow = dataRow - min(dataRow);
    dataRow = dataRow/max(dataRow);
%     dataRow = normalize(dataRow);
    dataRow = dataRow*0.8 - i;
    plotData(i,startIndex:endIndex) = dataRow;
end

plot(plotData');
