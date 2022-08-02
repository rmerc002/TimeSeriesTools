function [n_max] = getPaaFactor(timeseries, subsequenceLength)

user_patience = 2; % user patience is set to 2 seconds
margin = 0.1;
low_th = user_patience - margin;
high_th = user_patience + margin;
while length(timeseries) < 10000
    timeseries = cat(1,timeseries, timeseries);
end
%n_max = min([10000, length(timeseries)]);
n_max = 10000;
initial = n_max;
step = 5000;

basetimeseries = timeseries(1:n_max); 
tic;
[~] = SimMat(basetimeseries, subsequenceLength);
original_time = toc;
ns = [];
while true
    if ((original_time >= low_th) && (original_time <= high_th)) || (length(ns) > 2 && n_max == ns(length(ns)-2))
        break;
    
    elseif original_time < low_th
        n_max = n_max + step;
        if length(timeseries) <= n_max
            %n_max = length(timeseries);
            %break;
            timeseries = cat(1,timeseries, timeseries);
        end
        basetimeseries = timeseries(1:n_max); 
        tic;
        [~] = SimMat(basetimeseries, subsequenceLength);
        original_time = toc;
        
    elseif original_time > high_th
        n_max = n_max - step;
        basetimeseries = timeseries(1:n_max); 
        tic;
        [~] = SimMat(basetimeseries, subsequenceLength);
        original_time = toc;
        
    end
    ns(end+1) = n_max;
    if n_max == initial
       step = ceil(step / 2);
    end
end
    
warning('Maximum length to be computed in user patience time is set to >> %d',n_max);   
warning('Time to compute the maximum length >> %d',original_time);
%paa = ceil(length(timeseries)/n_max);
%warning('Downsampling rate is set to %d',paa);
