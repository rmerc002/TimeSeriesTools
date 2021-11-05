function [ts,dts, Seq, pChoice,names] = DNAFilesToTimeseries(path)
%%% path to a directory with DNA text files. path DOES NOT end with '/'

%%%

path = strcat(path,'/');
files = dir(strcat(path,'*.txt'));
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
    temp = [sum(txt=='a'), sum(txt=='c'), sum(txt=='g'), sum(txt=='t')];
    charSums(index,:) = temp;
end

[charSumss,ind] = sort(sum(charSums,1),'descend');
disp(charSumss)
disp(ind)
pChoice = zeros(1,4);

pChoice(1,ind(1)) = 1;
pChoice(1,ind(2)) = -charSumss(1)*1/charSumss(2);
pChoice(1,ind(3)) = 2;
pChoice(1,ind(4)) = -charSumss(3)*2/charSumss(4);
disp(pChoice);

ts = {};
dts = {};
for index = 1:size(files,1)
    dts_temp = zeros(1,length(Seq{index}));
    dts_temp(Seq{index} == 'a') = pChoice(1);
    dts_temp(Seq{index} == 'c') = pChoice(2);
    dts_temp(Seq{index} == 'g') = pChoice(3);
    dts_temp(Seq{index} == 't') = pChoice(4);
    dts_temp(dts_temp == 0) = randn(1,sum(dts_temp == 0))*0.00001;

    dts{index} = dts_temp;
    
    ts_temp = zeros(1,length(dts_temp));
    for index2=2:length(dts_temp)+1
        ts_temp(index2) = ts_temp(index2-1) + dts_temp(index2-1);
    end
    ts{index} = ts_temp(1:end-1);
    
end
end

