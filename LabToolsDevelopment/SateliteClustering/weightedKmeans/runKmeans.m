function [preds, centroids, kmeansOutput, clustdist] = runKmeans(data, row, col, nclust, seg)
    
    %Inputs : 
    %               1) data : Input time series data
    %               2) row_count : Output image row count
    %               3) col_count : Output image col count
    %               4) nclust : Number of clusters (Optional)
    %               5) seg : Segment size (Optional)
    % By default : 
    %               1) number_of_clusters  = 5;
    %               2) segment_size  = 1;

    clc;
    rng(0);
    data = zscore(data,0,1);
    
    if (nargin < 4)  nclust = 5; end;
    if (nargin < 5)  seg = 1; end;
    
    kmeansOutput = nan(row, col);
    colors = [1 0 0;... % red for 0
              0 1 0;... % green for 1
              0 0 1;... % blue for 2 (0 0 1)
              0 0 0; ... %black for 3
              0.4940 0.1840 0.5560; ... %purple for 4
              %1 0 1; % ... magneta for 5
              %0 1 1; % ... cyan for 6
              ];
    
    
    % No splitting and no ANR applied 
    [preds, centroids, ~, clustdist] = kmeansfold(data, nclust);
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
    %figure(1)
    %imagesc(kmeansOutput)
    %colormap(colors);
    %set(gca,'DataAspectRatio',[1 1 1])
    
    
    
    
    
    
    %%===========================================================
    %%===========================================================
    %{
    %nclust = 5;
    apply_anr = false;
    splitting = false;
    split_factor = 0;
    %seg = 10 ;
    %loc_x = readmatrix('/Users/maryam/Documents/phd/vert_ucr/results/ls_3_30_2022/ls_US_R001_whiteRoof/timeseries/lsma/lsma_band_2_sw_pixels_segment_xloc.xlsx');
    %loc_y = readmatrix('/Users/maryam/Documents/phd/vert_ucr/results/ls_3_30_2022/ls_US_R001_whiteRoof/timeseries/lsma/lsma_band_2_sw_pixels_segment_yloc.xlsx');
    %loc_x = readmatrix('/Users/maryam/Documents/phd/vert_ucr/results/wv_features/timeseries/US_R001_0002_MAT_argmax_singleband/US_R001_0002_MAT_argmax_singleband_band_0_sw_ixels_segment_xloc.xlsx');
    %loc_y = readmatrix('/Users/maryam/Documents/phd/vert_ucr/results/wv_features/timeseries/US_R001_0002_MAT_argmax_singleband/US_R001_0002_MAT_argmax_singleband_band_0_sw_pixels_segment_yloc.xlsx');
    %loc = cat(2,loc_x, loc_y);
    %data = readmatrix('/Users/maryam/Documents/phd/vert_ucr/results/ls_3_30_2022/US_R001_whiteRoof/lsma_band_2_sw_timeseries.csv');
    %data = readmatrix('/Users/maryam/Documents/phd/vert_ucr/results/wv_features/timeseries/US_R001_0002_MAT_argmax_singleband/US_R001_0002_MAT_argmax_singleband_band_0_sw_timeseries.xlsx');
    %data = readmatrix('/Users/maryam/Documents/phd/vert_ucr/results/wv_features/timeseries/US_R001_0029_MAT_argmax_singleband/US_R001_0029_MAT_argmax_singleband_band_0_sw_timeseries.csv');
    %data = readmatrix('/Users/maryam/Documents/phd/vert_ucr/results/ls_transient_events/timeseries/burning_swaped/burning_swaped_band_7_sw_timeseries.csv');
    %main_data = data;
    %data = cat(1,main_data,rand(2310720,size(main_data,2)));
    %data = zscore(data,0,1);
    %init_centroids = readmatrix('/Users/maryam/Documents/phd/vert_ucr/results/wv_features/timeseries/US_R001_0002_MAT_argmax_singleband/sw_kmeansmask_band_0_nc_5.csv');
    %data = data_spatial_corn_wheat(:,1:end-2);
    %data = data_spatial_corn_wheat(:,end-1:end);

    %data = data_spatial_USR0040012(:, 1:end-2);
    %data = zscore(data,0,1);
    %[totrow, totcol] = size(data);
    %for ii=1:size(data,1)
    %    data(ii,:) = zscore(data(ii,:));
    %end
    %%Korea lsma
    %row = 288;
    %col = 269;
    %%lsma white roof
    %row = 251;
    %col = 317;
    %%single band US-R001-0002
    %row = 82;
    %col = 132;
    %%single band US-R001-0029
    %row = 181;
    %col = 159;
    %%burning man
    %row=58;
    %col=59;
    %%SITS_Fold1
    %row = 100;
    %col = 10;
    %%US_R004_0012
    %row = 88;
    %col = 136;

    if splitting
        %data_anr = data_up;
        row_up = floor(row/split_factor);
        row_dn = floor(row/split_factor);
        kmeansOutput_up = nan(row_up, col);
        kmeansOutput_dn = nan(row_dn, col);
    else
        data_anr = data;
        kmeansOutput = nan(row, col);
    end

    if apply_anr
        [data_isClusterable] = ANR(data_anr);
        if splitting
            data_filtered = data_up(data_isClusterable==1,:);
        else
            data_filtered = data(data_isClusterable==1,:);
        end
    end

    if splitting
        if apply_anr
            [preds_up_filtered, centroids_up] = kmeansfold(data_filtered, nclust);
            [~,preds_dn] = pdist2(centroids_up([1,3:end],:),data_dn,'euclidean','Smallest',1);
            preds_up = nan(size(data_up,1),1);
            ic = 1;
            for ip = 1:size(data_up,1)
                if data_isClusterable(ip)
                    preds_up(ip) = preds_up_filtered(ic);
                    ic = ic + 1;
                else
                    preds_up(ip) = nclust;
                end
            end
        else
            [preds_up, centroids_up] = kmeansfold(data_up, nclust);
            [~,preds_dn] = pdist2(centroids_up([1,3:end],:),data_dn,'euclidean','Smallest',1);
        end
        shift = 1; % Sliding Window
        ts_i = 1;
        ii = 1;
        while ii < row_up - seg - 1
            jj =1;
            while jj < col - seg - 1
                kmeansOutput_up(ii:ii+shift, jj:jj+shift) = preds_up(ts_i);
                ts_i = ts_i + 1;
                jj = jj +shift;
            end
            ii = ii + shift;
        end

        shift = 1; % Sliding Window
        ts_i = 1;
        ii = 1;
        while ii < row_dn - seg - 1
            jj =1;
            while jj < col - seg - 1
                kmeansOutput_dn(ii:ii+shift, jj:jj+shift) = preds_dn(ts_i);
                ts_i = ts_i + 1;
                jj = jj +shift;
            end
            ii = ii + shift;
        end

    else
        if apply_anr
            [preds_filtered, centroids,sumd, clustdist] = kmeansfold(data_filtered, nclust);
            preds = nan(size(data,1),1);
            ic = 1;
            for ip = 1:size(data,1)
                if data_isClusterable(ip)
                    preds(ip) = preds_filtered(ic);
                    ic = ic + 1;
                else
                    preds(ip) = nclust;
                end
            end
        else
            [preds, centroids, sumd, clustdist] = kmeansfold(data, nclust);
            %train_centroids = nan(size(centroids_US_R001_0002,1), size(data,2));
            %size(train_centroids)
            %for nn=1:size(centroids_US_R001_0002,1)
            %    NewCt = linspace(1, size(centroids_US_R001_0002,2), size(data,2));
            %    train_centroids(nn,:) = interp1(linspace(1, size(centroids_US_R001_0002,2), size(centroids_US_R001_0002,2)), centroids_US_R001_0002(nn,:), NewCt); 
            %end
            %[~,preds] = pdist2(train_centroids,data,'euclidean','Smallest',1);

        end
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
    end





    if apply_anr
        colors = [1 0 0;... % red for 0
                  0 1 0;... % green for 1
                  0 0 1;... % blue for 2
                  0 0 0; ... %black for 3
                  0.4940 0.1840 0.5560; ... %purple for 4
                  %1 0 1; % ... magneta for 5
                  %0 1 1; % ... cyan for 6
                  1 1 1;    % white for 5 (Don't know)
                  ];
    else
        colors = [1 0 0;... % red for 0
                  0 1 0;... % green for 1
                  0 0 1;... % blue for 2 (0 0 1)
                  0 0 0; ... %black for 3
                  0.4940 0.1840 0.5560; ... %purple for 4
                  %1 0 1; % ... magneta for 5
                  %0 1 1; % ... cyan for 6
                  ];
    end
    colors_dn = [1 0 0;... % red for 0
      0 1 0;... % green for 1
      0 0 1;... % blue for 2
      0 0 0; ... %black for 3
      0.4940 0.1840 0.5560; ... %purple for 4
      ];

    if splitting
        figure(1)
        subplot(2,1,1)
        imagesc(kmeansOutput_up)
        colormap(colors);
        figure(2)
        subplot(2,1,2)
        imagesc(kmeansOutput_dn)
        colormap(colors_dn);
    else
        figure(1)
        imagesc(kmeansOutput)
        colormap(colors);
        set(gca,'DataAspectRatio',[1 1 1])
    end
    %}
    %%================================================
    %%================================================
    %%================================================
    %{
    figure(2)
    for i=1:nclust
        subplot(nclust,1,i);
        plot(centroids_up(i,:),'Color',colors(i,:));
    end
    %}

    function [preds, cent, sumd, clustdist] = kmeansfold(data, nc)
    [preds,cent, sumd, clustdist] = kmeans(data,nc);
    %[preds,cent] = kmedoids(data,nc,'Distance',@mixdistance);
    %[preds,cent,sumd, clustdist] = kmedoids(data,nc);
    preds = preds-1;

    end

end