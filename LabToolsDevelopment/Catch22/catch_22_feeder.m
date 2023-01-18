%% Call The function to produce the figures.
% This function takes a UCR format dataset and returns the Catch 22 features
% of that dataset.
% The UCR format is [label, time-series].
%                   [label, time-series]. 
%                           .
%                           . 
function val=catch_22_feeder(data)
        [rows,~]=size(data);
        vals=[];
        for i =1:rows
            label=data(i,1);
            x=data(i,2:end)';
            [new,~]=catch22_all(x);
            l=[label,new'];
            vals=[vals;l];    
        end
        val=vals;
end