function weightedkmeans(data, row, col, weight, nclust, seg)
    clc;
    %Inputs : 
    %               1) data : Input time series data
    %               2) row_count : Output image row count
    %               3) col_count : Output image col count
    %               4) nclust : Number of clusters (Optional)
    %               5) seg : Segment size (Optional)
    % By default : 
    %               1) number_of_clusters  = 5;
    %               2) segment_size  = 1;

    if (nargin < 5)  nclust = 5; end;
    if (nargin < 6)  seg = 1; end;
    
    
    data_p1 = data(:,1:end-2);
    data_p2 = data(:,end-1:end);
    
    [~, ~, ~, clustdist_part1] = runKmeans(data_p1, row, col, nclust, seg);
    [~, ~, ~, clustdist_part2] = runKmeans(data_p2, row, col, nclust, seg);
    
    
    %US_ls_whiteroof
    %row = 251;
    %col = 317;
    %weight = 2;
    %seg = 10;
    %test_corn_wheat
    %row = 100;
    %col = 10;
    %weight = 3;
    
    sumdistances = clustdist_part1 + weight * clustdist_part2;
    [~,preds] = min(sumdistances,[],2);
    preds = preds - 1;

    kmeansOutput = nan(row, col);

    colors = [1 0 0;... % red for 0
              0 1 0;... % green for 1
              0 0 1;... % blue for 2 (0 0 1)
              0 0 0; ... %black for 3
              0.4940 0.1840 0.5560; ... %purple for 4
              ];      

    count = 0;      
    shift = 1; % Sliding Window
    ts_i = 1;
    ii = 1;
    while ii < row - seg
        jj =1;
        while jj < col - seg
            kmeansOutput(ii:ii+shift, jj:jj+shift) = preds(ts_i);
            ts_i = ts_i + 1;
            jj = jj +shift;
        end
        ii = ii + shift;
    end

    f = figure(1);
    imagesc(kmeansOutput)
    colormap(colors);
    set(gca,'DataAspectRatio',[1 1 1])
    %saveas(f,sprintf('WV_US_R004_0012_%d.svg',weight))
    
end