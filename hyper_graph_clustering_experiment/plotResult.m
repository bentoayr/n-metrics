function newFig = plotResult(legendOff,cate_config,ave,std,algSet,xtag,ytag,displayOrder)
if nargin ==7
    displayOrder = 1:length(algSet.algEnable);
end
plotSet.lineWidth = 2.5;
plotSet.markerSize = 10;
plotSet.tickFontSize = 15;
plotSet.legendFontSize = 24;
plotSet.labelFontSize = 20;
plotSet.frameWidth = 1.5;
paran = length(cate_config);
algn = length(algSet.algEnable);
   newFig = figure('NumberTitle', 'off','Name', [xtag, ' ', ytag],'color','white');
    valideLegendSet = [];
    for ii = 1:algn
        algk = displayOrder(ii);
        if algSet.algEnable(algk)==0,continue;end
            valideLegendSet = [valideLegendSet,algk];
            plot(cate_config,ave(1:paran,algk),'LineWidth', plotSet.lineWidth, ...
                    'Color', algSet.algColor{algk}, ...
                    'LineStyle', algSet.algLineStyle{algk}, ...
                    'Marker', algSet.algMarker{algk}, ...
                    'MarkerSize', plotSet.markerSize...
                    );
            hold on;
            dx = abs(diff(cate_config));
            if numel(dx)==0, dx=1;end
            wid2 = dx(1) / 2;
            devWid = 0.15;

            hold on;
    end
    Xmin = min(cate_config); Xmax = max(cate_config);

    Ymin = min(min(ave(1:paran,valideLegendSet)))-0.01; 
    Ymax = max(max(ave(1:paran,valideLegendSet)))+0.01;
    x1 = Xmin-wid2*devWid;x2 = Xmax+wid2*devWid;
    if x1==x2, x1=double(x1)-0.05; x2 = double(x2)+0.05; end
    y1 = min(0.9,Ymin);y2 = (Ymax-max(0,Ymin))*0.05+Ymax;
    if y1==y2, y1=y1-0.05; y2 = y2+0.05; end
    if cate_config(end)<=cate_config(1)
        set(gca,'xdir','reverse');
    end
    
    axis([x1 x2  y1 y2]);
    set(gca,'FontName','Arial','FontSize',plotSet.tickFontSize);
    if legendOff==0
        hLegend = legend(algSet.algNameSetDisplay{valideLegendSet},'Location','Best');
    end
    axis square;
    hold on;

    xlabel(xtag,'Fontsize',plotSet.labelFontSize,'FontName','Arial','FontWeight','bold');hold on;
    ylabel(ytag,'Fontsize',plotSet.labelFontSize,'FontName','Arial','FontWeight','bold');hold on;
    if cate_config(end)<=cate_config(1)
        set(gca,'xtick',cate_config(end:-1:1));
    else
        set(gca,'xtick',cate_config);
    end
    if y1<1 && y2<=1.5
        ytickNum = floor(y1*10)/10:0.1:y2;
        for yk = 1:length(ytickNum)
            ytickText{yk} = num2str(ytickNum(yk),'%.1f');
        end
        set(gca,'ytick',ytickNum);
        set(gca,'YTickLabel',ytickText);
    end
    set(gca,'linewidth',plotSet.frameWidth);
    box on;