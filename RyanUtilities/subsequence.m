classdef subsequence

    properties
        subsequenceWithContext;
        contextLength;
        subsequenceLength;
        length;
        tsMin;
        tsMax;
        tsStartIndexContext;
        tsStartIndexSubsequence;
    end
    methods
        function obj = subsequence(ts, subsequenceIndex, subsequenceLength, contextLength)
            obj.subsequenceLength = subsequenceLength;
            obj.contextLength = contextLength;
            obj.subsequenceWithContext = obj.extractSubsequenceFromTS(ts, subsequenceIndex);
            obj.length = 2*contextLength + subsequenceLength;
            obj.tsMin = min(ts);
            obj.tsMax = max(ts);

            obj.tsStartIndexContext = subsequenceIndex - contextLength;%%%Can be negative
            obj.tsStartIndexSubsequence = subsequenceIndex;
        end
    
        function ssi = subsequenceStartIndex(obj)
            ssi = obj.contextLength+1;
        end

        function ssi = subsequenceEndIndex(obj)
            ssi = obj.contextLength+obj.subsequenceLength;
        end
    
        function ss = subsequenceWithoutContext(obj)
            startIndex = obj.subsequenceStartIndex;
            endIndex = startIndex + obj.subsequenceLength-1;
            ss = obj.subsequenceWithContext(startIndex:endIndex);
        end
    
        function subsequenceWithContext = extractSubsequenceFromTS(obj, ts, subsequenceStartIndex)
            subsequenceWithContext = nan(2*obj.contextLength+obj.subsequenceLength,1);
            
            startIndex = subsequenceStartIndex - obj.contextLength;
            startSSIndex = 1 + max(0, 1+(-1*startIndex));
            startIndex = max(1, startIndex);
    
            endIndex = subsequenceStartIndex + obj.subsequenceLength + obj.contextLength-1;
            endSSIndex = length(subsequenceWithContext) + min(0, length(ts)-endIndex);
            endIndex = min(length(ts), endIndex);
    
            subsequenceWithContext(startSSIndex:endSSIndex) = ts(startIndex:endIndex);
        end
    end
end
