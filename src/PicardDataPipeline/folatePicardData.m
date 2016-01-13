configPicard;
v2struct();
disp('running folatePicardData');
outputDir1 = [outputDir filesep 'folatePicardData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

hugoFolate = {'MTR','BHMT2','MAT1A','MAT2A','MAT2B','GNMT','AHCY','AHCYL1',...
    'AHCYL2','CHDH','BAD','CBS','CTH','GCLC','GCLM','GSS','SHMT1','SHMT2',...
    'SARDH','GLDC','ALDH1L1','ALDH1L2','MTHFD1','MTHFS','GART','ATIC','TYMS',...
    'MTHFR','DHFR'};

selectExpPicardData(hugoFolate, configStruct, outputDir1);