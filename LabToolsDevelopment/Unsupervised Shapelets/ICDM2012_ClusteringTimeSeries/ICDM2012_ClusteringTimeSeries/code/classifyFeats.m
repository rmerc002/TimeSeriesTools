function accuracy = classifyFeats(feats,labels)
    n = size(feats,1);
    correct = 0;    
    for i = 1:n
        bsf = inf;
        minI = -1;
        for j = 1:n
            if i == j
                continue;
            end
            d = ED(feats(i,:),feats(j,:));
            if d < bsf
                bsf = d;
                minI = j;
            end
        end
        if(minI==-1)
            mn=0.02;
        end
        if labels(minI) == labels(i)
            correct = correct + 1;
        end
    end

    accuracy = correct/n;

end