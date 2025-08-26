function f = subplotRSOM(FigureNumSinglePlots,FigureNumSubplot, opt)
    arguments
        FigureNumSinglePlots
        FigureNumSubplot
        opt.FastScanningAxis          {mustBeMember(opt.FastScanningAxis, {'x','y'})}       = 'y';
        opt.Visible            char   {mustBeMember(opt.Visible, {'on', 'off'})}            = 'on';             % if off, figures become invisible        
        opt.FigurePosition                                                                  = [];
    end
    
    opt.defaultFigVis = get(groot, 'DefaultFigureVisible');
    set(groot, 'DefaultFigureVisible', 'off');
    
    if ishandle(FigureNumSubplot)
        close(FigureNumSubplot)
    end
    switch  opt.FastScanningAxis
        case 'x'
            set(0, 'CurrentFigure', FigureNumSinglePlots(1));
            ax = gca;
            view(90,90);
            ax.XAxis.Direction = "reverse";
            FastScanningAxis = 3;
            SlowScanningAxis = 2;
        case 'y'
            set(0, 'CurrentFigure', FigureNumSinglePlots(1));
            ax = gca;
            view(0,90);
            ax.XAxis.Direction = "normal";
            FastScanningAxis = 2;
            SlowScanningAxis = 3;
    end
    
    j = 0;
    for k = FigureNumSinglePlots
        j = j+1;
        set(0, 'CurrentFigure', k);
        ax(j) = gca;
    end
    f = figure(FigureNumSubplot);
    for j = 1:3 
        ax_copy(j) = copyobj(ax(j),f);
        switch j
            case 1
                subplot(3,3,[1,2,3],ax_copy(j))
            case FastScanningAxis
                subplot(3,3,[4,5,7,8],ax_copy(j))
            case SlowScanningAxis
                subplot(3,3,[6,9],ax_copy(j))
        end    
    end
    f.Visible = opt.Visible;
    if strcmp(opt.defaultFigVis, 'on')
        f = figure(FigureNumSubplot);
    end
    if ~isempty(opt.FigurePosition)
        f.Position(5-numel(opt.FigurePosition):end) = opt.FigurePosition;
    end
    set(groot, 'DefaultFigureVisible', opt.defaultFigVis);
end


