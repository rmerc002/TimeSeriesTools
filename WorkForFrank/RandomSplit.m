%%
% Returns a binary tree traversal of the integers from 1 to n
function [IDX] = RandomSplit(Range)
IDX = [1];
indices = 2:length(Range);

IDX = [IDX, indices(randperm(length(indices)))];

end
