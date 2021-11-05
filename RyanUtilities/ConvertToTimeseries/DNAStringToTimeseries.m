function [ts,dts, pChoice] = DNAStringToTimeseries(dnaString)
%%% path to a directory with DNA text files. path DOES NOT end with '/'

%%%


txt = strtrim(dnaString);
txt = lower(txt);
Seq = txt;
temp = [sum(txt=='a'), sum(txt=='c'), sum(txt=='g'), sum(txt=='t')];
charSums = temp;

[charSumss,ind] = sort(sum(charSums,1),'descend');
disp(charSumss)
disp(ind)
pChoice = zeros(1,4);

pChoice(1,ind(1)) = 1;
pChoice(1,ind(2)) = -charSumss(1)*1/charSumss(2);
pChoice(1,ind(3)) = 2;
pChoice(1,ind(4)) = -charSumss(3)*2/charSumss(4);
disp(pChoice);


    dts = zeros(1,length(Seq));
    dts(Seq == 'a') = pChoice(1);
    dts(Seq == 'c') = pChoice(2);
    dts(Seq == 'g') = pChoice(3);
    dts(Seq == 't') = pChoice(4);
    dts(dts == 0) = randn(1,sum(dts == 0))*0.00001;
    
    ts = zeros(1,length(dts));
    for index2=2:length(dts)+1
        ts(index2) = ts(index2-1) + dts(index2-1);
    end

end

