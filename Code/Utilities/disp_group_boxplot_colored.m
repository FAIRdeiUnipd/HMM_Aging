function disp_group_boxplot_colored(data,group)

% data: n x m x d matrix
% boxplot will be grouped into m groups, each of 6 boxplots describing the
% distribution of the d values 

temp = [1:length(group)];
temp2 = repmat(temp',1,size(data,1))+[0 0.1 0.2 0.3 0.4 0.5];
positions = reshape(temp2',size(data,1)*length(group),1);
boxplot(data','positions', positions',...
    'OutlierSize',7,'Symbol','k+')
color = [0.5 0.8 0.3; 0.85 0.20 0.098; 0.1010 0.4450 0.933; 0.343 0 0; 0.643 0 0; 0 0 1];

color = repmat(color,9,1);
h = findobj(gca,'Tag','Box');
for j = 1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.9);
end
c = get(gca, 'Children');
legend
hleg1 = legend(c(1:6),'S1','S2','S3','S4','S5','S6','Orientation','horizontal','NumColumns',1,'Location','eastoutside');
set(gca,'XTick',[]) ,set(gcf,'Position',get(0,'Screensize'))
set(gca,'FontSize',16,'FontWeight','bold')
xticklabels('')

