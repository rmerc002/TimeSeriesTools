function [weekends, weekdays] = PHAllLengths(profileIndices, subLenSeries)
    weekends = zeros(1,length(subLenSeries));
    weekdays = zeros(1,length(subLenSeries));
    timeIndices = 1:size(profileIndices,2);
    for lengthIndex = 1:length(subLenSeries)
       nonNan = ~isnan(profileIndices(lengthIndex,:));
       
       profileIndices(lengthIndex,~nonNan) = timeIndices(~nonNan);%quick approx by removing nan
       
       persevEnergy = zeros(1,size(profileIndices,2));
       persevEnergy(nonNan) = abs(timeIndices(nonNan) - profileIndices(lengthIndex,nonNan)) < 30*24;
       
       [weekdayTotal, weekdaySum, weekendTotal, weekendSum] = calculateItalianPowerFraction(profileIndices(lengthIndex,:), persevEnergy);
       
       %part of the quick fix
       weekdayTotal = weekdayTotal - sum(~nonNan);
       weekdaySum = weekdaySum - sum(~nonNan);
       weekendTotal = weekendTotal - sum(~nonNan);
       weekendSum = weekendSum - sum(~nonNan);
       
       weekends(lengthIndex) = weekendSum/weekendTotal;
       weekdays(lengthIndex) = weekdaySum/weekdayTotal;
    end
end