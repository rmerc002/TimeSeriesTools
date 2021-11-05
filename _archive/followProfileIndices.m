function [profileIndices2] = followProfileIndices(profileIndices)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
profileIndices2 = profileIndices;

%set best neighbors equal to smaller
%the smaller one will point to itself
for row=1:size(profileIndices,1)
    for col=1:size(profileIndices,2)
        col2 = profileIndices(row,col);
        if col2 == 0
            continue;
        end
%         if col == profileIndices(row,col2)
%             profileIndices2(row,col) = col;
%             profileIndices2(row,col2) = col;
%         end
        profileIndices2(row,col) = rabbitHole(col, profileIndices(row,:));
    end
    row
end

end

function index2 = rabbitHole(index, profileIndices1D)
    indexNN = profileIndices1D(index);
    if index == profileIndices1D(indexNN)
        index2 = min(index, indexNN);
    else
        index2 = rabbitHole(indexNN, profileIndices1D);
    end
end