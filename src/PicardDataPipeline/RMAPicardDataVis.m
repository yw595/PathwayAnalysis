configPicard;
v2struct();
disp('running RMAPicardDataVis');
outputDir1 = [outputDir filesep 'RMAPicardData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

load([outputDir1 filesep 'RMAPicardDataPermute.mat']);

tempNames1 = {'refNormScores','refMets','refMetNames','refPVals', 'permuteScoresMatrix'};
tempNames2 = {'T','Corr','Anova'};

for i=1:1%length(tempNames2)
    for k=1:length(tempNames1)
            eval([tempNames1{k} ' = ' tempNames1{k} tempNames2{i} ';']);
            eval(['titleString = RMAPermute' tempNames2{i} ';'])
    end
    for j=1:1%2
        permuteScoresMatrix = permuteScoresMatrix(:,:,j);
        if j==1
            titleString = [titleString 'permuteAll'];
        else
            titleString = [titleString 'permuteSig'];
        end
         
        xlabels = refMets;
        ylabelString = 'Norm Score';
        xvals = []; yvals = [];
        outputDir = outputDir1;
        for k=1:size(permuteScoresMatrix,1)
            disp(k)
            xvals = [xvals; k*ones(size(permuteScoresMatrix,2),1)];
            yvals = [yvals; permuteScoresMatrix(k,:)'];
        end
        
        figure('Visible','off');
        scatter(xvals,yvals,[],'b','MarkerFaceColor','b','SizeData',10);
        hold on;
        scatter(1:length(refNormScores),refNormScores,[],'r','MarkerFaceColor','r','SizeData',60);
        hold on;
        ylabel(ylabelString,'FontSize',20);
        ylim([1.3*min(yvals(:)) 1.3*max(yvals(:))]);
        title(titleString,'FontSize',20);
        set(gca,'XTickLabel',{});
        set(gca,'FontSize',20);
        line([0 length(xlabels)], [0 0]);
        saveas(gcf,[outputDir filesep titleString '.png']);
        close(gcf);
    end
end