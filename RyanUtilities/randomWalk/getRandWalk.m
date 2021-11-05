function [rwalk] = getRandWalk(length)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
steps = randi([-1,1],10,length);
rwalk = zeros(10,length);
for i = 2:length
    rwalk(:,i) = rwalk(:,i-1) + steps(:,i);
end
rwalk(2:end,:) = rwalk(2:end,:) - rwalk(2:end,end);

splicePointFound = false;
for i=2:499
    for i2=2:10
        if rwalk(1,i) == rwalk(i2,i)
            rwalk(1,i:end) = rwalk(i2,i:end);
            splicePointFound = true;
            break
        end
    end
    if splicePointFound == true
        break
    end
end
rwalk = rwalk(1,:);
rwalk = rwalk./max(abs(rwalk));
    
end

