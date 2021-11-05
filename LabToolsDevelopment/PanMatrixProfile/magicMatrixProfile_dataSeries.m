function  [magicMP, profileIndices,yAxis] = magicMatrixProfile_dataSeries(dataSeries, subLength)
if size(dataSeries,1) > size(dataSeries,2)
    dataSeries = dataSeries';
end

yAxis = 1:size(dataSeries,1);

magicMP = zeros(size(dataSeries)); %preallocate
profileIndices = zeros(size(dataSeries)); %preallocate
for dataIndex = 1:size(dataSeries,1)
    disp(dataIndex);
    
    data = normalize(dataSeries(dataIndex,:));
    
    [matrixProfile,profileIndex] = mpx(data',ceil(subLength/2),subLength);

    magicMP(dataIndex,1:length(matrixProfile)) = matrixProfile;
    profileIndices(dataIndex,1:length(profileIndex)) = profileIndex;
    
end

%%% remove values that are zero at the end of the matrix profile due to
%%% subLenth
    magicMP(:,end-(subLength):end) = nan;

end