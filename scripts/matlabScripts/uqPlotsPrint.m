function [ ] = uqPlotsPrint(location, fileName, precisionX, precisionY)
    %% Ticks format
    if nargin <= 2
      precisionX = 0;
      precisionY = 0;
    end
    
    ax = gca;
    ax.XRuler.TickLabelFormat = ['%.' int2str(precisionX) 'f'];
    ax.YRuler.TickLabelFormat = ['%.' int2str(precisionY) 'f'];

    %% TIKZ
%     matlab2tikz([location '/tex/' fileName '.tex'],'showInfo',false,...
%                   'height','*\figureheight','width','\figurewidth',...
%                   'parseStrings', true);
%     matlab2tikz([location '/tex/' fileName '.tex'],'showInfo',false,...
%                   'height','5cm','width','6cm',...
%                   'parseStrings', false);
    
    %% PDF (cropped)
%     tmpVar = [location '/pdf/' fileName '.pdf'];
%     print(gcf,'-dpdf',tmpVar);
%     [ status, cmdout ] = system(sprintf("pdfcrop %s %s",tmpVar,tmpVar));
    
    %% Tranparent Background
    set(gca, 'Color', 'none');

    %% PNG
    figname = [location '/png/' fileName];
%     export_fig(figname,'-transparent','-r300','-png');   
    print(gcf,[location '/png/' fileName '.png'],'-dpng','-r300');
    
end
