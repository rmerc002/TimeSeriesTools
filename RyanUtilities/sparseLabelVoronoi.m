function voronoiLabels = sparseLabelVoronoi(sparseLabels)
    indices = 1:length(sparseLabels);
    nonZeroIndices = indices(sparseLabels>0);

    voronoiLabels = zeros(length(sparseLabels),1);
    voronoiLabels(1:nonZeroIndices(1)) = sparseLabels(nonZeroIndices(1));

    for ii = 1:length(nonZeroIndices)-1
        midIndex = ceil((nonZeroIndices(ii+1)+nonZeroIndices(ii))/2);
        voronoiLabels(nonZeroIndices(ii):midIndex) = sparseLabels(nonZeroIndices(ii));
        voronoiLabels(midIndex:nonZeroIndices(ii+1)) = sparseLabels(nonZeroIndices(ii+1));
    end
    voronoiLabels(nonZeroIndices(end):length(sparseLabels)) = sparseLabels(nonZeroIndices(end));
end