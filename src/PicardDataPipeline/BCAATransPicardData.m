configPicard;
v2struct();
disp('running BCAATransPicardData');
outputDir1 = [outputDir filesep 'BCAATransPicardData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

hugoBCAATrans = {'SLC38A1','SAT1','SLC7A5','SLC7A8','SLC7A7','SLC7A6','SLC3A2','SLC7A5',...
    'SLC43A1','SLC43A2','SLC3A2'};

selectExpPicardData(hugoBCAATrans, configStruct, outputDir1);