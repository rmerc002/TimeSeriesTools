function [rwalk] = getRandWalkBasic(length)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
steps = randn(1,length);
rwalk = zeros(1,length);
for i = 2:length
    rwalk(:,i) = rwalk(:,i-1) + steps(:,i);
end
rwalk(2:end,:) = rwalk(2:end,:) - rwalk(2:end,end);

    
end

