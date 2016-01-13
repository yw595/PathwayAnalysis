configPicard;
v2struct();
disp('running adHoc');
outputDir1 = [outputDir filesep 'adHoc'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

firstHeteroplasmies = zeros(length(PicardHeteroplasmies),1);
for j=1:length(uniqueHeteroplasmies)
    for k=1:length(PicardHeteroplasmies)
        if strcmp(PicardHeteroplasmies{k},uniqueHeteroplasmies{j})
            firstHeteroplasmies(k)=1;
            break;
        end
    end
end
firstHeteroplasmies = find(firstHeteroplasmies);

load([outputDir filesep 'fluxPicardData' filesep 'fluxResults.mat']);

outputPrefixes = {'FALCON','EFlux','iMAT','GIMME'};

% for i=1:length(outputPrefixes)
%     for k=1:2
%         for j=1:length(firstHeteroplasmies)-1
%             titleString = outputPrefixes{i};
%             if k==1
%                 titleString = [titleString ' RPMIMed'];
%             else
%                 titleString = [titleString ' ShlomiMed'];
%             end
%             titleString = [titleString ' ' uniqueHeteroplasmies{j} ' Vs ' uniqueHeteroplasmies{j+1}];
%             fluxDist1 = fluxResults{i,firstHeteroplasmies(j),k};
%             fluxDist2 = fluxResults{i,firstHeteroplasmies(j+1),k};
%             [~, sortIdxs] = sort(abs(fluxDist1-fluxDist2)); sortIdxs = sortIdxs(end:-1:end-9);
%             fluxDist1 = fluxDist1(sortIdxs); fluxDist2 = fluxDist2(sortIdxs); fluxDist = [fluxDist1 fluxDist2];
%             figure('Visible','off'); b=bar(1:10,fluxDist); %b(1).FaceColor='red'; b(2).FaceColor='blue';
%             ylabel('Flux','FontSize',20);
%             xlabels = origRecon2.rxns(sortIdxs);
%             ylim([-1.3*max(max(abs(fluxDist))) 1.3*max(max(abs(fluxDist)))]);
%             for l=1:length(xlabels)
%                 text(l,-1.3*max(max(abs(fluxDist))),xlabels{l},'Rotation',90,'FontSize',10);
%             end
%             set(b(1),'FaceColor','b');
%             set(b(2),'FaceColor','r');
%             legend(b,{uniqueHeteroplasmies{j},uniqueHeteroplasmies{j+1}});
%             title(titleString);
%             disp([outputDir1 filesep strrep(titleString,' ','') '.png']);
%             saveas(gcf,[outputDir1 filesep strrep(titleString,' ','') '.png']);
%             close(gcf);
%         end
%     end
% end
% 
refFluxResults = fluxResults;
load([outputDir filesep 'fluxPicardData' filesep 'permuteFluxResultsSafe2.mat']);

for i=1:1
    for k=1:2
        titleString = [outputPrefixes{i} 'Safe'];
        if k==1
            titleString = [titleString ' RPMIMedPermute'];
        else
            titleString = [titleString ' ShlomiMedPermute'];
        end
        
        xlabels = {};
        ylabelString = 'Flux Diff';
        xvals = []; yvals = [];
        outputDir = outputDir1;
        refDiffDistDisplay = [];
        for j=1:length(firstHeteroplasmies)-1
            refFluxDist1 = refFluxResults{i,firstHeteroplasmies(j),k};
            refFluxDist2 = refFluxResults{i,firstHeteroplasmies(j+1),k};
            [~, sortIdxs] = sort(abs(refFluxDist1-refFluxDist2)); sortIdxs = sortIdxs(end:-1:1);
            refDiffDist = refFluxDist1-refFluxDist2; refDiffDist = refDiffDist(sortIdxs);
            refDiffDistDisplay = [refDiffDistDisplay; refDiffDist];
            for l=1:length(permuteFluxResults)
                disp([num2str(i) ' ' num2str(k) ' ' num2str(j) ' ' num2str(l)])
                fluxResults = permuteFluxResults{l};
                fluxDist1 = fluxResults{i,firstHeteroplasmies(j),k};
                fluxDist2 = fluxResults{i,firstHeteroplasmies(j+1),k};
                %[~, sortIdxs] = sort(abs(fluxDist1-fluxDist2)); sortIdxs = sortIdxs(end:-1:1);
                fluxDist1 = fluxDist1(sortIdxs); fluxDist2 = fluxDist2(sortIdxs); diffDist = fluxDist1-fluxDist2;
                
                xvals = [xvals; j*(1:length(diffDist))'];
                yvals = [yvals; diffDist];
            end
        end
        figure('Visible','off');
        scatter(xvals,yvals,[],'b','MarkerFaceColor','b','SizeData',10);
        hold on;
        scatter(1:length(refDiffDistDisplay),refDiffDistDisplay,[],'r','MarkerFaceColor','r','SizeData',60);
        hold on;
        ylabel(ylabelString,'FontSize',20);
        ylim([1.3*min(yvals(:)) 1.3*max(yvals(:))]);
        title(titleString,'FontSize',20);
%         for k=1:length(xlabels)
%             disp(k)
%             text(k,-.1*max(abs(yvals(:))),xlabels{k},'Rotation',90,'FontSize',20);
%         end
        set(gca,'XTickLabel',{});
        set(gca,'FontSize',20);
        for j=1:length(firstHeteroplasmies)-1
            line([j*length(diffDist) j*length(diffDist)], [1.3*min(yvals(:)) 1.3*max(yvals(:))]);
        end
        saveas(gcf,[outputDir filesep titleString '.png']);
        close(gcf);
    end
end