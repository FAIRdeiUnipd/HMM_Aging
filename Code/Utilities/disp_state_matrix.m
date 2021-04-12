function disp_state_matrix(M,state_list,map_range)

K = length(state_list);
imagesc(M)
hold on
for j=0:(K+1)
    plot([0 K+1] - 0.5,[j j] + 0.5,'k','LineWidth',1)
    plot([j j] + 0.5,[0 K+1] - 0.5,'k','LineWidth',1)
end
caxis(map_range)
colorbar
xticks(1:K)
yticks(1:K)
set(gca,'fontsize',14,'fontweight','bold')
xlabel('To state','Fontsize',14,'FontWeight','bold')
ylabel('From state','Fontsize',14,'FontWeight','bold')
set(gca,'TickLength',[0,0])
set(gcf,'Position',get(0,'Screensize'));
hold off

end