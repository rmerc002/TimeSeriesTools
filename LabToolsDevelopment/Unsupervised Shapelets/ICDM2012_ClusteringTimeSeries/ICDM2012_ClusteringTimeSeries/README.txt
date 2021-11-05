-> The "Trace" dataset used in the paper is included in the sampleDataset folder
-> All the codes are included in code folder

----------Instruction to run the code----------

-> Run 'featureSelect.m'
   it will generate,
     1. acc1,
     2. loc1,
     3. ftr1,
     4. the extracted shapelets are saved in the folder "acc1"
     5. the 'shapelets' are marked with red in the time series
     
-> Run 'distCorr.m'
  it will generate the plot to show the change in Rand index,
     1. with the addition of the 'shapelets' (green curve)
     2. Rand index using ground truth (blue curve)
     3. Rand index when entire time series is used (red line)
 
-----------------------------------------------.


-> The two synthetic datasets, rock dataset and eloctronic usage dataset used
in the paper are not publicly available.