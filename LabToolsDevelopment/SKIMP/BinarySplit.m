%%
% Returns a binary tree traversal of the integers from 1 to n
function [IDX] = BinarySplit(Range)
n = length(Range);
% Preallocate memory to store the traversed indices and the remaining intervals
IDX       = [];
Intervals = {};

IDX(1) = 1; % We always begin by explore the first integer
Intervals(1) = {[2,n]}; % After exploring the first integer, we begin splitting the interval 2:n

counter = 0;
while(~isempty(Intervals))
    fprintf("counter: %d\n",counter);
    counter = counter + 1;
    
    lb = Intervals{1}(1);
    ub = Intervals{1}(2);
    mid = floor((lb + ub) / 2);
    Intervals(1) = [];
    
    IDX(end+1) = mid;
    
    if(lb == ub)
        continue;
    else
        [L,R] = split(lb,ub,mid);
        if(~isempty(L))
            Intervals{end+1} = L;
        end
        if(~isempty(R))
            Intervals{end+1} = R;
        end
    end
end


function [L,R] = split(lb,ub, m)
    if(lb == m)
        L = [];
        R = [m+1, ub];
    elseif ub == m
        L = [lb m-1];
        R = [];
    else
        L = [lb m-1];
        R = [m+1 ub];
    end
end
            
            
end
