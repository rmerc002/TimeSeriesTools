function [featureAnnotations] = loadMouseFeatureAnnotation(fileName)
    fid = fopen(fileName);
    featureAnnotations = {};
    for ii = 1:4
        line = string(fgets(fid));
        featureAnnotations{ii} = split(line, ",");
    end
    fclose(fid);
end