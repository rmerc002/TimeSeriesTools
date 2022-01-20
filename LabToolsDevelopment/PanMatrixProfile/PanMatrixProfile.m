function PanMatrixProfile(tsA, startLength, endLength, numLengths)
    [magicMP, profileIndices,subLenSeries] = magicMatrixProfile(tsA, startLength, endLength, numLengths);
    visualizeMMPAB(tsA, magicMP, profileIndices, [], [], [], subLenSeries, subLenSeries, "subLength", "");
end