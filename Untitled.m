%get maps of mouseID to LineNum, QLength, Sex and Tissue
[junk1, junk2, rawMetadata] = xlsread( ...
    'C:\\Users\\Yiping Wang\\My Documents\\Lab Work\\HD\\6_month_allelic_series_data_mRNA_with_metadata2014_21.xlsx','mRNA Metadata');
mouseIDToInfo = containers.Map;
for i=2:size(rawMetadata,1)
    aMouseID = rawMetadata{i,1};
    aLineNum = i;
    aQLength = str2num(rawMetadata{i,4}(2:end));
    aSex = rawMetadata{i,6};
    aTissue = rawMetadata{i,3};
    infoStruct = struct('lineNum',aLineNum,'QLength',aQLength,'sex',aSex,'tissue',aTissue);
    mouseIDToInfo(aMouseID) = infoStruct;
end

[junk1, ~, rawData] = xlsread( ...
    'C:\\Users\\Yiping Wang\\My Documents\\Lab Work\\HD\\6_month_allelic_series_data_mRNA_with_metadata2014_21.xlsx','mRNA FPKM data');
ensemblIDToInfo = containers.Map;
FCortexMask = makeOrderedMask(mouseIDToInfo, mouseIDList,'sex','F','tissue','cortex');
FLiverMask = makeOrderedMask(mouseIDToInfo, mouseIDList,'sex','F','tissue','Liver');
FStriatumMask = makeOrderedMask(mouseIDToInfo, mouseIDList,'sex','F','tissue','striatum');
MCortexMask = makeOrderedMask(mouseIDToInfo, mouseIDList,'sex','M','tissue','cortex');
MLiverMask = makeOrderedMask(mouseIDToInfo, mouseIDList,'sex','M','tissue','Liver');
MStriatumMask = makeOrderedMask(mouseIDToInfo, mouseIDList,'sex','M','tissue','striatum');
for i=2:size(rawData,1)
    i
    ensemblID = rawData{i,1};
    expressionVals1 = [];expressionVals2 = [];expressionVals3 = [];
    expressionVals4 = [];expressionVals5 = [];expressionVals6 = [];
    QVals1 = [];QVals2 = [];QVals3 = [];QVals4 = [];QVals5 = [];QVals6 = [];
    for j=2:size(rawData,2)
        expressionVal = rawData(i,j);
        mouseID = rawData{1,j};       
        if(strcmp(mouseIDToInfo(mouseID).sex,'F'))
            if(strcmp(mouseIDToInfo(mouseID).tissue,'cortex'))
                expressionVals1(end+1) = expressionVal{1};QVals1(end+1) = mouseIDToQLength(mouseID);
            elseif(strcmp(mouseIDToTissue(mouseID),'Liver'))
                expressionVals2(end+1) = expressionVal{1};QVals2(end+1) = mouseIDToQLength(mouseID);
            elseif(strcmp(mouseIDToTissue(mouseID),'striatum'))
                expressionVals3(end+1) = expressionVal{1};QVals3(end+1) = mouseIDToQLength(mouseID);
            end
        elseif(strcmp(mouseIDToSex(mouseID),'M'))
            if(strcmp(mouseIDToTissue(mouseID),'cortex'))
                expressionVals4(end+1) = expressionVal{1};QVals4(end+1) = mouseIDToQLength(mouseID);
            elseif(strcmp(mouseIDToTissue(mouseID),'Liver'))
                expressionVals5(end+1) = expressionVal{1};QVals5(end+1) = mouseIDToQLength(mouseID);
            elseif(strcmp(mouseIDToTissue(mouseID),'striatum'))
                expressionVals6(end+1) = expressionVal{1};QVals6(end+1) = mouseIDToQLength(mouseID);
            end
        end
    end
    [rho pval] = corr();
    [rho1 pval1] = corr(QVals1', expressionVals1', 'type', 'Spearman'); pvals1corr(end+1) = pval1; rhos1(end+1) = rho1;
    [rho2 pval2] = corr(QVals2', expressionVals2', 'type', 'Spearman'); pvals2corr(end+1) = pval2; rhos2(end+1) = rho2;
    [rho3 pval3] = corr(QVals3', expressionVals3', 'type', 'Spearman'); pvals3corr(end+1) = pval3; rhos3(end+1) = rho3;
    [rho4 pval4] = corr(QVals4', expressionVals4', 'type', 'Spearman'); pvals4corr(end+1) = pval4; rhos4(end+1) = rho4;
    [rho5 pval5] = corr(QVals5', expressionVals5', 'type', 'Spearman'); pvals5corr(end+1) = pval5; rhos5(end+1) = rho5;
    [rho6 pval6] = corr(QVals6', expressionVals6', 'type', 'Spearman'); pvals6corr(end+1) = pval6; rhos6(end+1) = rho6;
    [h pval1] = ttest2(expressionVals1(QVals1==20)', expressionVals1(QVals1==175)'); pvals1t(end+1) = pval1;
    [h pval2] = ttest2(expressionVals2(QVals2==20)', expressionVals2(QVals2==175)'); pvals2t(end+1) = pval2;
    [h pval3] = ttest2(expressionVals3(QVals3==20)', expressionVals3(QVals3==175)'); pvals3t(end+1) = pval3;
    [h pval4] = ttest2(expressionVals4(QVals4==20)', expressionVals4(QVals4==175)'); pvals4t(end+1) = pval4;
    [h pval5] = ttest2(expressionVals5(QVals5==20)', expressionVals5(QVals5==175)'); pvals5t(end+1) = pval5;
    [h pval6] = ttest2(expressionVals6(QVals6==20)', expressionVals6(QVals6==175)'); pvals6t(end+1) = pval6;
end

save('readRNASeqDataMitoMiner.mat');