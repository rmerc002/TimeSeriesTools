clc;
clear;
warning off;


%%%%rand index for algo
D = 'C:\ftrSlcn\ICDM2012_ClusteringTimeSeries\sampleDataset';
data = load('C:\ftrSlcn\ICDM2012_ClusteringTimeSeries\sampleDataset\Trace.txt');
cls = data(:,1);

nr = size(data,1);
clsPrev = ones(nr,1);
clstr=4;



RDS=[];
for ftrNum = 1:6,  %;
    DIS=[];
    fid = fopen('C:\ftrSlcn\ICDM2012_ClusteringTimeSeries\sampleDataset\ftr1.txt','r');
    for i=1:ftrNum
        ftr = fgetl(fid);
        r=regexp(ftr,'\s','split');
        FTR = [];
        l=length(r)-1;
        for j=1:l,
            FTR(j)=str2double(r{j});
        end
        dist = computeMatrix(FTR,data);
        dist = dist/sqrt(l);
        DIS = [DIS dist];
    end
    fclose(fid);



    %%%%    d/(r1+r2)
    bst=inf;
    DIS = real(DIS);

    for i=1:20
        [IDX,C,SUMD]=kmeans(DIS,clstr,'EmptyAction','singleton');
        d=pdist(C);
        rds = [];   %%radius of each circle
        

        if(sum(SUMD)<bst)
            bst=sum(SUMD);
            clsCurr = IDX;
            RI = RandIndex(cls,IDX); %% rand index using ground truth
            RI_=RandIndex(clsPrev,clsCurr); %% rand index with the addition of shapelet
        end
    end
    clsPrev = clsCurr;
    RDS(ftrNum)=bst;
    RDSRI(ftrNum)=RI;
    RDSRI_(ftrNum)=RI_;
end

%%%%rand index using entire time series %%%

bst=inf;
for i=1:20
    [IDX,C,SUMD]=kmeans(data(:,2:end),clstr);
    d=pdist(C);
    rds = [];   %%radius of each circle
    
    if(sum(SUMD)<bst)
        bst=sum(SUMD);
        RIts = RandIndex(cls,IDX);
    end
end

figure;
plot(RDSRI,'b-*');
hold on;
plot(1-RDSRI_,'g-*');
plot([1 6],[RIts RIts],'r-');

