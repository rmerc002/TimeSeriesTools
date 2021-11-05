function extractTextGivenIndicesAndSubLength(text1, text2,downsampleStep, index1, index2, subLength)
    subsequence1 = text1(index1*downsampleStep-(downsampleStep-1):index1*downsampleStep-(downsampleStep-1)+subLength*downsampleStep-(downsampleStep-1));
    subsequence2 = text2(index2*downsampleStep-(downsampleStep-1):index2*downsampleStep-(downsampleStep-1)+subLength*downsampleStep-(downsampleStep-1));
    disp(subsequence1);
    disp(subsequence2);
end