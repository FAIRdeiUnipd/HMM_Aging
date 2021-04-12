function disp_mean_std_area_colored(avg_data,std_data,thr,label_x,isZscore)

if isZscore == 0
    
    color = [0 1 1; 1 0 1; 1 1 0; 0 0 1; 0 1 0; 1 0 0];
    
    h1=figure
    for t = 1:length(thr)
        s_thr = thr(t);
        
        subplot(3,3,t)
        MEAN_VALUE      = squeeze(avg_data(:,:,t));
        MEAN_PLUS_1SD   = squeeze(avg_data(:,:,t))+squeeze(std_data(:,:,t));
        MEAN_MINUS_1SD  = squeeze(avg_data(:,:,t))-squeeze(std_data(:,:,t));
        
        values = [1:length(MEAN_VALUE)];
        
        for k = 1: size(MEAN_VALUE,1)
            hold on
            h=fill([values fliplr(values)],[MEAN_MINUS_1SD(k,:) fliplr(MEAN_PLUS_1SD(k,:))],color(k,:),'linestyle','none');
            set(h,'facealpha',.1);
            plot(values,MEAN_VALUE(k,:),'o-','Color',color(k,:),'MarkerFaceColor',color(k,:));
        end
        
        hold off
        set(gca,'fontsize',10,'fontweight','bold')
        xlim([1 length(label_x)])
        xticks([1:length(label_x)])
        xticklabels({label_x{:,1}})
        c = get(gca, 'Children');
        %         hleg1 = legend(c(1:2:end),'S1','S2','S3','S4','S5','S6','Orientation','horizontal','NumColumns',3,'Location','South');
        hleg1 = legend(c(end-1:-2:1),'S1','S2','S3','S4','S5','S6','Orientation','horizontal','NumColumns',3,'Location','South');
        %     legend boxoff
        %     legend('mean \pm 1sd','mean')
        title(['sparsity thr=' num2str(s_thr)])
        set(gcf,'Position',get(0,'Screensize'))
        grid on
    end
    
elseif isZscore == 1
    
    color = [0 1 1; 1 0 1; 1 1 0; 0 0 1; 0 1 0; 1 0 0];
    
    h1=figure
    for t = 1:length(thr)
        s_thr = thr(t);
        subplot(3,3,t)
        hold on
        MEAN_VALUE      = squeeze(avg_data(:,:,t));
        for k = 1:size(MEAN_VALUE,1)
            plot(MEAN_VALUE(k,:),'o-','Color',color(k,:),'MarkerFaceColor',color(k,:));
        end
        set(gca,'fontsize',10,'fontweight','bold')
        xlim([1 length(label_x)])
        xticks([1:length(label_x)])
        xticklabels({label_x{:,1}})
        c = get(gca, 'Children');
        legend('S1','S2','S3','S4','S5','S6','Orientation','horizontal','NumColumns',3,'Location','South');
        %     legend boxoff
        %         plot(1.5*ones(1,11),'k-','LineWidth',1.5)
        %         plot(-1.5*ones(1,11),'k-','LineWidth',1.5)
        hold off
        title(['sparsity thr=' num2str(s_thr)])
        set(gcf,'Position',get(0,'Screensize'))
        grid on
    end
    
elseif isZscore == 2
    
    color = [0 1 1; 1 0 1; 1 1 0; 0 0 1; 0 1 0; 1 0 0; 0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0; 0 0 0.5; 0 0.5 0;...
        0.5 0 0; 0.2 0.2 0.2; 1 0.2 0.2; 0.2 0.2 1];
    
    h1=figure
    hold on
    MEAN_VALUE      = avg_data;
    for k = 1:size(MEAN_VALUE,1)
        plot(MEAN_VALUE(k,:),'o-','Color',color(k,:),'MarkerFaceColor',color(k,:));
    end
    set(gca,'fontsize',10,'fontweight','bold')
    xlim([1 length(label_x)])
    xticks([1:length(label_x)])
    xticklabels({label_x{:,1}})
    c = get(gca, 'Children');
    legend('d12','d13','d14','d15','d16','d23','d24','d25','d26','d34','d35','d36','d45','d46','d56'...
        ,'Orientation','horizontal','NumColumns',5,'Location','South');
    %     legend boxoff
    %         plot(1.5*ones(1,11),'k-','LineWidth',1.5)
    %         plot(-1.5*ones(1,11),'k-','LineWidth',1.5)
    hold off
    set(gcf,'Position',get(0,'Screensize'))
    grid on
    
elseif isZscore ==3
    
    color = [0 1 1; 1 0 1; 1 1 0; 0 0 1; 0 1 0; 1 0 0];
    
    h1=figure  
    
    MEAN_VALUE      = avg_data;
    MEAN_PLUS_1SD   = avg_data+std_data;
    MEAN_MINUS_1SD  = avg_data-std_data;
    
    values = [1:length(MEAN_VALUE)];
    
    for k = 1: size(MEAN_VALUE,1)
        hold on
        h=fill([values fliplr(values)],[MEAN_MINUS_1SD(k,:) fliplr(MEAN_PLUS_1SD(k,:))],color(k,:),'linestyle','none');
        set(h,'facealpha',.1);
        plot(values,MEAN_VALUE(k,:),'o-','Color',color(k,:),'MarkerFaceColor',color(k,:),'MarkerSize',10);
    end
    
    hold off
    set(gca,'fontsize',14,'fontweight','bold')
    xlim([1 length(label_x)])
    xticks([1:length(label_x)])
    xticklabels({label_x{:,1}})
    c = get(gca, 'Children');
    %         hleg1 = legend(c(1:2:end),'S1','S2','S3','S4','S5','S6','Orientation','horizontal','NumColumns',3,'Location','South');
    hleg1 = legend(c(end-1:-2:1),'S1','S2','S3','S4','S5','S6','Orientation','horizontal','NumColumns',3,'Location','South');
    %     legend boxoff
    %     legend('mean \pm 1sd','mean')
    set(gcf,'Position',get(0,'Screensize'))
    grid on
    
elseif isZscore ==4
    
%     color = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.494 0.184 0.556; 0.466 0.674 0.188; 0.3010 0.7450 0.933];
%    color = [0 0.447 0.741; 0.543 0 0; 0.85 0.325 0.098; 0.3010 0.7450 0.933; 0.929 0.694 0.125; 0.466 0.674 0.188];
    color = [0 0 1; 0.643 0 0; 0.343 0 0; 0.1010 0.4450 0.933; 0.85 0.20 0.098; 0.5 0.8 0.3];
   
    MEAN_VALUE      = avg_data;
    values = [1:length(MEAN_VALUE)];
    
    for k = 1: size(MEAN_VALUE,1)
        hold on
        plot(values,MEAN_VALUE(k,:),'o-','Color',color(k,:),'MarkerFaceColor',color(k,:),'MarkerSize',9,'LineWidth',1.3);
    end
    
    hold off
    set(gca,'fontsize',16,'fontweight','bold')
     ylim([-2 2])
    xlim([1 length(label_x)])
    xticks([1:length(label_x)])
    xticklabels({label_x{:,1}})
    xtickangle(45)
    c = get(gca, 'Children');
    legend('S1','S2','S3','S4','S5','S6','Orientation','horizontal','NumColumns',1,'Location','eastoutside');
    set(gcf,'Position',get(0,'Screensize'))
end

end
