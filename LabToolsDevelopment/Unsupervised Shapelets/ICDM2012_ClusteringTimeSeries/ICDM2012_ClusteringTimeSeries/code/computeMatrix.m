function dis = computeMatrix(feats, data)
n = size(data,1);
F = size(feats,1);
dis(1:n,1:F) = 0;
for j = 1:n
    for i = 1:F
        [loc(j,i) dis(j,i)] = findNN(data(j,2:end),feats(i,:));
    end
end