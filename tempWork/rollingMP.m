function [mp] = rollingMP(data, w)
%ROLLINGMP Summary of this function goes here
%   Detailed explanation goes here
data = data;
mp = inf(length(data),length(data));
maxVal = sqrt(2*w);
for row = 1:length(data)-w+1
    for col = row:length(data)-w+1
        if abs(row-col) <= ceil(w/2)
           mp(row,col) = maxVal;
           mp(col,row) = maxVal;
        else
            distance = inf;
            for rollOffset = 1:w
                tempDistance = norm(zscore(data(row:row+w-1)) - zscore(circshift(data(col:col+w-1),rollOffset)));
                distance = min(distance, tempDistance);
                
            end
            if distance < mp(row,col) 
                  mp(row,col) = distance;
                  mp(col,row) = distance;
            end
        end
    end
end

%normalize for length
mp = mp./(sqrt(2*w));
mp(end-w+2:end,:) = [];
mp(:,end-w+2:end) = [];

mp(eye(length(data)-w+1)==1) = 0;

end

