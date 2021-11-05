function [ts, dts, Seq, charValMapping, names] = TextToTimeseries(path,charValMapping)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% path = '/Users/cdslug/Dropbox/Education/UCR/Research/0_Projects/Rastamat/Paper/GraphFiles/DNAWalkComparison/';
if path(end-3:end) == '.txt' 
    files = dir(path);
else
    path = strcat(path,'/');
    files = dir(strcat(path,'*.txt'));
end
chars = '';
charSums = zeros(size(files,1),4);
Seq = {};
names = {};
for index = 1:size(files,1)
    file = files(index);
%     file.name;
    disp(file.name);
    names{index} = file.name;
    txt = fileread(strcat(path,file.name));
    txt = strtrim(txt);
    txt = lower(txt);
    Seq{index} = txt;
    chars = unique(strcat(chars,txt));
end

if isempty(charValMapping)
    charSums = zeros(1,length(chars));
    for fileIndex = 1:size(files,1)
        for charIndex = 1:length(chars)
            charSums(charIndex) = sum(Seq{fileIndex} == chars(charIndex));
        end 
    end

    [charSums,sortedIndices] = sort(charSums,'descend');
    chars = chars(sortedIndices);
    pChoice = zeros(1,length(chars));

    for i=1:2:length(chars)
        pChoice(i) = (i+1)/2;

        if i == length(chars)
            pChoice(i) = (i+1)/2;%last character if chars is odd
            break;
        end
        pChoice(i+1) = -pChoice(i)*charSums(i)/charSums(i+1);
    end

    charValMapping = containers.Map();
    for index = 1:length(chars)
        charValMapping(chars(index)) = pChoice(index);
    end
end


%%% Use the Choice %%%
ts{size(files,1)} = [];
dts{size(files,1)} = [];

for fileIndex = 1:size(files,1)
    dts{fileIndex} = zeros(1,length(Seq{fileIndex}));
    charValMappingKeys = charValMapping.keys();
    for index=1:length(charValMappingKeys())
        indexChar = charValMappingKeys(index);
        indexChar = indexChar{1};
        dts{fileIndex}(Seq{fileIndex} == indexChar) = charValMapping(indexChar);
    end

    ts{fileIndex} = cumtrapz(dts{fileIndex});
end

