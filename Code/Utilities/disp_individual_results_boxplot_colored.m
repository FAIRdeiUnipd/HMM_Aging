function disp_individual_results_boxplot_colored(data,Nyoung,widths,label_xtick,label_x)

% data: n x m x d matrix
% boxplot will be grouped into m groups, each of 6 boxplots describing the
% distribution of the d values 

K       = size(data,2);
temp    = [1:K];
temp2   = repmat(temp',1,2)+[0 0.2]; 
positions = reshape(temp2',2*K,1);

boxplot([data(1:Nyoung,1),data(Nyoung+1:end,1),data(1:Nyoung,2),data(Nyoung+1:end,2),...
    data(1:Nyoung,3),data(Nyoung+1:end,3),data(1:Nyoung,4),data(Nyoung+1:end,4),...
    data(1:Nyoung,5),data(Nyoung+1:end,5),data(1:Nyoung,6),data(Nyoung+1:end,6)],'positions', positions',...
    'Widths',widths,'OutlierSize',7,'Symbol','k+')
color = [1 0 0; 0 0 1]; %% first color is the last to be assigned to the boxplot

color = repmat(color,K,1);
h = findobj(gca,'Tag','Box');
for j = 1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5);
end
c = get(gca, 'Children');
legend
hleg1 = legend(c(1:2),'YOUNG','OLD','Orientation','horizontal','NumColumns',3,'Location','North');
legend boxoff, set(gca,'TickLength',[0,0]), set(gcf,'Position',get(0,'Screensize'))
set(gca,'FontSize',15,'FontWeight','bold')
xticks([1.1:1:6.1])
xticklabels(label_xtick)
xlabel(label_x)

