function [mp] = naiveMP(data, w)
%ROLLINGMP Summary of this function goes here
%   Detailed explanation goes here
mp = inf(1,length(data));

for row = 1:length(data)-w+1
%     if mod(row,100) == 0
%        disp(row); 
%     end
    for col = 1:length(data)-w+1
        if abs(col-row) > ceil(w/2)
            distance = norm(zscore(data(row:row+w-1)) - zscore(data(col:col+w-1)));
            mp(row) = min(mp(row),distance);
        end

        
    end
end

%normalize for length
mp = mp./(sqrt(2*w));
mp(end-w+2:end) = nan;

end

