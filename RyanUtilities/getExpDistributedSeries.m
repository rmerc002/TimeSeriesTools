function subLenSeries = getExpDistributedSeries(startLen, endLen,numSteps)
%%% Purpose: Reduce space. matrix profile distances are pretty similar from
%%% from one subLen to the next. By trial and error, adding the square root
%%% of the current subLen seems to produce a good distribution.
%     subLenSeries = [];
%     index = startLen;
%     
%     while index < endLen
%         subLenSeries = [subLenSeries,index];
%         index = index + max(1,ceil(sqrt(index)));
%     end
%     
%     subLenSeries = [subLenSeries,endLen];

powerMin = log10(startLen);
powerMax = log10(endLen);
% numSteps = 80;
powerStep = (powerMax-powerMin)/numSteps;
powers = powerMin:powerStep:powerMax;
subLenSeries = unique(ceil(power(10,powers)));

end