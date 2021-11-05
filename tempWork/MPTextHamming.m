function mp = MPTextHamming(A,B,subLength)

mp = inf(1,length(A) - subLength + 1);
exclusionWidth = 0;
if strcmp(A,B) == true
    exclusionWidth = subLength;
end

for iA = 1:length(A) - subLength + 1
   for iB = 1:length(B) - subLength + 1 
       if abs(iA - iB) >= exclusionWidth
           tempDist = pdist([A(iA:iA+subLength-1);B(iB:iB+subLength-1)],'hamming')*subLength;
           mp(iA) = min(mp(iA), tempDist);
       end
   end
end

end