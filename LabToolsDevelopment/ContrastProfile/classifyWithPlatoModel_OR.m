function [ypred] = classifyWithPlatoModel_OR(X, platoModel)
    ypred = -1*ones(size(X,1),1);
    numClasses = size(platoModel,1);

    for sampleIndex = 1:size(X,1)
        sample = X(sampleIndex,:);
        closestClassIndex = NaN;
        minDist = inf;
        %%%test each plato of each class
        for posClassIndex = 1:numClasses
            classRMSE = 0;
            for negClassIndex = 1:numClasses
                if posClassIndex == negClassIndex
                    continue;
                end
                
                
                platoStruct = platoModel{posClassIndex, negClassIndex};
                numPlatos = platoStruct.bestK;
                minDist2 = inf;
                for platoIndex = 1:numPlatos
                    plato = platoStruct.platos(platoIndex,:);
                    distanceProfile = real(MASS_V2(sample,plato));
                    dist = min(distanceProfile)/sqrt(2*length(plato));
                    dist = min(1, dist);
                    if dist < minDist2
                        minDist2 = dist;
                    end
                end
                classRMSE = classRMSE + minDist2^2;
            end
           classRMSE = sqrt(classRMSE/numPlatos);
           if classRMSE < minDist
               minDist = classRMSE;
               closestClassIndex = posClassIndex;
           end
        end
        ypred(sampleIndex) = closestClassIndex; 
    end
end