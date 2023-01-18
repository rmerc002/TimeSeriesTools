function subsequences = generateRandWalkShapeSet(shapeLength, numShapes, noiseAmp, baseShapeLength)
baseShape = getRandWalk(baseShapeLength);
shapeHiRes = interp1(1:baseShapeLength, baseShape, linspace(1,baseShapeLength,shapeLength));

subsequences = repmat(shapeHiRes,numShapes,1);
subsequences = subsequences + randn(numShapes, shapeLength)*noiseAmp;
end