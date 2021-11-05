function [weekdayTotal, weekdaySum, weekendTotal, weekendSum] = calculateItalianPowerFraction(profileIndices, persevEnergy)
weekdayTotal = 0;
weekdaySum = 0;
weekendTotal = 0; 
weekendSum = 0;
for hourIndex = 1:24*7:size(profileIndices,2)
    %Sunday
    if size(profileIndices,2) - hourIndex  >= 24
        weekendTotal = weekendTotal + 24;
        weekendSum = weekendSum + sum(persevEnergy(hourIndex:hourIndex + 24 - 1));
    else
%         disp('Sunday');
        remainingHours = size(profileIndices,2) - hourIndex;
        weekendTotal = weekendTotal + remainingHours;
        weekendSum = weekendSum + sum(persevEnergy(hourIndex: hourIndex + remainingHours - 1));
        break;
    end
    
    %Monday to Friday
    if size(profileIndices,2) - hourIndex >= 24*6
        weekdayTotal = weekdayTotal + 24*5;
        weekdaySum = weekdaySum + sum(persevEnergy(hourIndex + 24: hourIndex + 24 + 5*24 - 1));
    else
%         disp('weekday')
        remainingHours = size(profileIndices,2) - hourIndex - 24;
        weekdayTotal = weekdayTotal + remainingHours;
        weekdaySum = weekdaySum + sum(persevEnergy(hourIndex + 24: hourIndex + 24 + remainingHours - 1));
        break
    end
    
    %Saturday
    if size(profileIndices,2) - hourIndex >= 24*7
        weekendTotal = weekendTotal + 24;
        weekendSum = weekendSum + sum(persevEnergy(hourIndex + 24*6: hourIndex + 24*6 + 24 - 1));
    else
%         disp('Saturday')
        weekendTotal = weekdayTotal + size(profileIndices,2) - hourIndex - 24*6;
        weekendSum = weekendSum + sum(persevEnergy(hourIndex + 24*6:hourIndex + 24*6 +remainingHours - 1));
        break
    end
end
end
    