function plotQueryIndexNN(ts, queryIndex, m, KNN)
    query = ts(queryIndex:queryIndex+m-1);
    plotQueryNN(ts, query, KNN);
end