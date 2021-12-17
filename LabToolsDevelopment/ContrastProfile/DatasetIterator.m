classdef DatasetIterator
    properties
        path
        datasetIndex
        datasets
        numDatasets
    end
    methods
        function obj = DatasetIterator(datasetPath)
            obj.path = datasetPath;
            obj.datasetIndex = 0;
            tempFiles = dir(obj.path);
            obj.datasets = tempFiles(3:end);
            obj.numDatasets = length(obj.datasets);
        end
        function [obj, datasetIndex] = nextDatasetIndex(obj, steps)
            obj.datasetIndex = obj.datasetIndex + steps;
            if obj.datasetIndex > obj.numDatasets || obj.datasetIndex < 1
                obj.datasetIndex = 0;
            end
            datasetIndex = obj.datasetIndex;
        end
        function datasetPath = getFilePath(obj, datasetIndex)
            if datasetIndex > 0 && datasetIndex <= obj.numDatasets
                datasetPath = obj.datasets(datasetIndex);
            else
                datasetPath = "";
            end 
        end
        function datasetName = getDatasetName(obj, datasetIndex)
            if datasetIndex > 0 && datasetIndex <= obj.numDatasets
                datasetPath = obj.datasets(datasetIndex);
            else
                datasetPath = "";
            end 
        end
        function [samples, labels] = loadData(obj, datasetIndex, splitType, extension)
            samples = [];
            labels = [];

            datasetPath = obj.datasets(datasetIndex).folder;
            datasetName = obj.datasets(datasetIndex).name;
            if datasetName(1) == '.'
                return;
            end
        
            fullFileName = fullfile(datasetPath, datasetName,sprintf("%s_%s%s",datasetName, splitType, extension));
            fprintf(1, 'Now reading %s\n', datasetName);

            if ~exist(fullFileName, 'file') 
                fprintf("\n%%%%%%%%%% ERROR: skipping: %s\n\n", datasetName);
               return; 
            end
        
            data = importdata(fullFileName);
        
            samples = data(:,2:end);
            labels = data(:,1);
        end
    end
end