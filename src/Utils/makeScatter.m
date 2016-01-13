function makeScatter(xvals, yvals, titleString, xlabels, ylabelString, outputDir)
    figure('Visible','off');
    scatter(xvals,yvals,[],'b','MarkerFaceColor','b','SizeData',40);
    ylabel(ylabelString,'FontSize',20);
    ylim([-.1*max(yvals(:)) 1.3*max(yvals(:))]);
    title(titleString,'FontSize',20);
    for j=1:length(xlabels)
        text(j,-.1*max(yvals(:)),xlabels{j},'Rotation',90,'FontSize',20);
    end
    set(gca,'XTickLabel',{});
    line([0 length(xlabels)], [0 0]);
    saveas(gcf,[outputDir filesep titleString '.png']);
    close(gcf);
end