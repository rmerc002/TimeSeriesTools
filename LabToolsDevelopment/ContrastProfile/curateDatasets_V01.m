function curateDatasets(curatedDataPath)
%%% filename = "";
%%% descriptor
%%% Train/Test
%%% type = "segmented";
%%% needsTranspose = false;
%%% headerlines = 0;
%%% classIndex = 1;
%%% sampleStartIndex = 2;
%%% sampleEndIndex = 1025;

            
    %%% V01: reads all dirs from UCR archive     
    ucrArchivePath = "D:\Datasets\UCRArchive_2018\UCRArchive_2018";
    dirinfo = dir(ucrArchivePath);
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories

    datasets = {};
    for K = 1 : length(dirinfo)
        thisDataSetDir = dirinfo(K).name;
        if thisDataSetDir(1) == '.' 
            continue
        end
        if thisDataSetDir(2) == '.'
            continue
        end
        %         datasets{end+1} = {"D:\Datasets\UCRArchive_2018\UCRArchive_2018\Mallat\Mallat_TRAIN.tsv","Mallat","train","segmented",false,0,1,2,1025};
%         datasets{end+1} = {"D:\Datasets\UCRArchive_2018\UCRArchive_2018\Mallat\Mallat_TEST.tsv","Mallat","test","segmented",false,0,1,2,1025};
        disp(ucrArchivePath);
        disp(thisDataSetDir);
        disp(thisDataSetDir + "_Train.tsv");
        trainPath = fullfile(ucrArchivePath, thisDataSetDir, thisDataSetDir + "_TRAIN.tsv");
        datasets{end+1} = {trainPath,thisDataSetDir,"train","segmented",false,0,1,2,-1};
        
        trainPath = fullfile(ucrArchivePath, thisDataSetDir, thisDataSetDir + "_TEST.tsv");
        datasets{end+1} = {trainPath,thisDataSetDir,"test","segmented",false,0,1,2,-1};
    end

    
    
    %%% Use a different function of each set of data
    %%%   -This allows for looking back to see how it was created
    %%%   -Do not modify after it has been finalized and used in a
    %%%   presentation.
    %%%   -I need to make a comment about not modifying, maybe set permissions.
    datasetCreation_V01(datasets, curatedDataPath);
    
    
end