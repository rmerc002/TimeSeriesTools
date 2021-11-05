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

            
            
    ucrArchivePath = "D:\Datasets\UCRArchive_2018\UCRArchive_2018";

    datasets = {};
    dirnames = {'Mallat','BirdChicken', 'BME','ChlorineConcentration','CricketX','Crop','ECG200','EOGHorizontalSignal','FreezerSmallTrain'};
%     dirnames = {'ChlorineConcentration','CricketX','Crop','ECG200','EOGHorizontalSignal','FreezerSmallTrain'};
%     dirnames = {'Mallat','Car','CinCECGTorso','ElectricDevices','FaceAll','Fish'};
    
%     dirinfo1 = dir(ucrArchivePath);
%     dirinfo1(~[dirinfo1.isdir]) = [];  %remove non-directories
%     for K = 1 : length(dirinfo1)
%         thisDataSetDir = dirinfo1(K).name;
    
    for K = 1 : length(dirnames)
        thisDataSetDir = dirnames{K};
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
    datasetCreation_V08(datasets, curatedDataPath);
    
    
end