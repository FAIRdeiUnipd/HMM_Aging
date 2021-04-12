function disp_mean_state_values(mean_state,IC_list,comp_network_names)

% mean_state = vector containing the mean activation of the state
colors = {[1,1,0,0.4],[0.9290,0.6940,0.1250,0.4],[1,0,0,0.4],...
    [0.6350 0.0780 0.1840,0.4],[0.4940 0.1840 0.5560,0.4],[1,0,1,0.4],...
    [0,1,1,0.4],[0 0.4470 0.7410,0.4],[0,0,1,0.4],...
    [0,1,0,0.4],[0.4660 0.6740 0.1880,0.4]};

nIC = length(IC_list);

% figure
axis([0.5 nIC+0.5 min(mean_state)+min(mean_state)/2 ...
    max(mean_state)+max(mean_state)/2])
hold on
pos_x = 0;
for nn = 1:1:length(comp_network_names)
    if nn == 1
        position = [pos_x  min(mean_state)+min(mean_state)/2 ...
            length(comp_network_names{nn,2})+0.5...
            abs(min(mean_state)+min(mean_state)/2)+max(mean_state)+max(mean_state)/2];
        rect= rectangle('Position',position);
        rect.FaceColor = colors{1,nn};
        pos_x = pos_x + length(comp_network_names{nn,2})+0.5;
    else
        position = [pos_x  min(mean_state)+min(mean_state)/2 ...
            length(comp_network_names{nn,2}) ...
            abs(min(mean_state)+min(mean_state)/2)+max(mean_state)+max(mean_state)/2];
        rect= rectangle('Position',position);
        rect.FaceColor = colors{1,nn};
        pos_x = pos_x + length(comp_network_names{nn,2});
    end
    
end
errorbar(mean_state,[],'.k','MarkerSize',30)
% plot(0.5:1:nIC+0.5,zeros(nIC+1),'r','LineWidth',2)
xlabel('Independent component ID')
set(gca,'FontSize',15)
xticks(1:1:length(IC_list))
xticklabels(IC_list)
xtickangle(90)
xlim([0.5 nIC+0.5]);
grid on
set(gca,'TickLength',[0,0])
set(gcf,'Position',get(0,'Screensize'));