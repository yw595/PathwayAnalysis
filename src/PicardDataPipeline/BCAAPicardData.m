configPicard;
v2struct();
disp('running BCAAPicardData');
outputDir1 = [outputDir filesep 'BCAAPicardData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

hugoBCAA = {'BCAT2','BCKDK','PPM1K','BCKDHA','BCKDHB','DLD','DBT','IVD','HADHA',...
    'ECHS1','MCCC1','MCCC2','AUH','OXCT1','OXCT2','HMGCL','ACAD8','HIBCH','HIBADH','ALDH6A1',...
    'ACADSB','HSD17B10','ACAT1','ACAT2','PCCB'};

pValsBinom = binomTest(hugoBCAA, configStruct);

selectExpPicardData(hugoBCAA, configStruct, outputDir1);