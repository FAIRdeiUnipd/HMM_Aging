function disp_basic_matrix(matrix,list_IC,map_range)

imagesc(matrix)
set(gca,'TickLength',[0 0])
set(gca,'XTick',[],'YTick',[])

caxis(map_range)
colormap jet
colorbar('Fontsize',14)

nIC = length(list_IC);

xticks(1:nIC)
xtickangle(90)
xticklabels(list_IC)
yticks(1:nIC)
yticklabels(list_IC)
% xlabel('Independent component ID','Fontsize',15,'FontWeight','bold')
% ylabel('Independent component ID','Fontsize',15,'FontWeight','bold')

set(gca,'FontSize',8)
set(gcf,'Position',get(0,'Screensize'));
% set(gca,'FontSize',8,'FontWeight','bold')
