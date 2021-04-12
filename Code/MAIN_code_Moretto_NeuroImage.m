% Analysis code to support the findings of
% "Unveil whole-brain dynamics in normal aging through Hidden Markov
% Models"
% by Moretto et al., Neuroimage, 2021
%
% Code written by Moretto M.

clc
clear all
close all
warning off

% Set path of folder with Code and Data
folder_path = fullfile(pwd,'..');

Code_path       = fullfile(folder_path,'Code');
Utilities_path  = fullfile(pwd,'Utilities');
Results_path    = fullfile(folder_path,'Results');
Figures_path    = fullfile(folder_path,'Figures');
HMM_res_path    = fullfile(Results_path,'HMM_6states_rep5_conf2');

addpath(genpath(fullfile(Utilities_path)))
addpath(genpath(fullfile(Code_path)))

%da togliere poi
addpath(genpath('/nfsd/nasfair/Users/Manuela/TOOLS/SOFTWARE/NIfTI_20140122'))
addpath(genpath('/nfsd/nasfair/Users/Manuela/TOOLS/SOFTWARE/HMM-MAR-master'))
addpath(genpath('/nfsd/nasfair/Users/Manuela/TOOLS/SOFTWARE/BCT/2019_03_03_BCT'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K           = 6;
set_par     = 2;
folder_state    = [num2str(K) 'states'];

% Load options for HMM
option_name     = ['option_conf_' num2str(set_par) '.mat'];
load(fullfile(HMM_res_path,option_name));

Vol         = size(X,2);    %no. volumes
N           = size(X,1);    % no. subjects
Nyoung      = length(name_subj_young_25_sel);
Nold        = length(name_subj_old_60);

% Load FE for all the repetitions
load(fullfile(HMM_res_path,'fe_rep.mat'),'table','ord_table','best_rep');
r_idx       = 1;
r           = best_rep(r_idx,1);
% Load HMM
hmm_name    = ['HMM_rep_' num2str(r) '_conf_' num2str(set_par) '.mat'];
load(fullfile(HMM_res_path,hmm_name));

% Name of the ICs
load(fullfile(Results_path,'names_IC_RSNs.mat'))
[comp_network_names,slctIC_origIC] = get_goodICs();
nIC = length(slctIC_origIC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AVERAGE LOG-LIKELIHOOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      AVERAGE LOG-LIKELIHOOD, AIC, BIC  ===================')

avLL = GammaavLL(hmm,Gamma,Xi,cell2mat(T));
disp(['average LL' num2str(avLL)])

n_params = K*nIC + K*nIC*(nIC+1)/2 + K*nIC*(nIC+1)/2 + 1 + K*K + K + K*size(Gamma,1);
disp(['Number of parameters ' num2str(n_params)])

AIC_LL = -2*avLL + 2* n_params;
disp(['AIC ' num2str(AIC_LL)])

BIC_LL = -2*avLL + n_params * log(N*Vol);
disp(['BIC ' num2str(BIC_LL)])
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 1 SUPPLEMENTARY MATERIAL: MEAN BOLD ACTIVATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      MEAN BOLD ACTIVATION  ===================')

for k = 1:K
    hold on
    disp_mean_values(hmm.state(k).W.Mu_W,lst_rsn_name,comp_network_names)
    plot(0.5:1:nIC+0.5,zeros(nIC+1),'r','LineWidth',2)
    set(gca,'fontsize',12,'fontweight','bold')
    title(['\fontsize{15} Mean BOLD activation, S' num2str(k)])
    exportgraphics(gcf,fullfile(Figures_path,'Mu',['Fig1_SI_S' num2str(k) '.png']))
    close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 2: MEAN BOLD ACTIVATION SPATIAL MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      MEAN BOLD ACTIVATION SPATIAL MAPS ===================')

tmp         = load_untouch_nii(fullfile(Results_path,['NEW_MASK_ZSCORE_ASSIGN_WM_CSF_MASKED.nii.gz']));
grp_mask    = double(tmp.img);
[rr,cc,zz]  = size(grp_mask);

for k = 1:K    
    disp(['State' num2str(k)])
    mean_act_mask = nan*zeros(size(grp_mask));
    mean_act_mask_zeros = zeros(size(grp_mask));
    for i = 1:rr
        for j = 1:cc
            for l = 1:zz
                idx_IC = grp_mask(i,j,l);
                if idx_IC==0
                    continue
                else
                    mean_act_mask(i,j,l) =  hmm.state(k).W.Mu_W(find(slctIC_origIC==idx_IC));
                    mean_act_mask_zeros(i,j,l) =  hmm.state(k).W.Mu_W(find(slctIC_origIC==idx_IC));
                end
            end
        end
    end
    
%     % 3D NIFTI FILE WITH ALL Mus OF 46 ICs TOGETHER
%     NII = create_4D_nii(fullfile(Results_path,'NEW_MASK_ZSCORE_ASSIGN_WM_CSF_MASKED.nii.gz'),mean_act_mask);
%     save_untouch_nii(NII,fullfile(Figures_path,'Mu',['MEAN_3D_K' num2str(k) '.nii.gz']));
%     clear NII
%     NII = create_4D_nii(fullfile(Results_path,'NEW_MASK_ZSCORE_ASSIGN_WM_CSF_MASKED.nii.gz'),mean_act_mask_zeros);
%     save_untouch_nii(NII,fullfile(Figures_path,'Mu',['MEAN_3D_K' num2str(k) '_zeros.nii.gz']));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CVs USING THE DIAGONAL OF THE PRECISION MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      CVs USING THE DIAGONAL OF THE PRECISION MATRIX ===================')

for k = 1:K
    disp(['State' num2str(k)])
    new_cv2(k,:)        = 100 * (sqrt(diag(hmm.state(k).Precision)).^(-1));
    mean_new_cv2(k,:)   = mean(new_cv2(k,:));
end
dist_cv_mean = new_cv2-mean_new_cv2;
save(fullfile(Results_path,'CVs','CV.mat'),'new_cv2','mean_new_cv2','dist_cv_mean')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 2 SUPPLEMENTARY MATERIAL: CVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
for k = 1:K
    subplot(3,2,k)
    hold on
    plot(ones(1,nIC)*mean_new_cv2(k,:),'-','LineWidth',2)
    stem(new_cv2(k,:),'r','filled','MarkerSize',7)
    hold off
    title(['CV S' num2str(k)])
    ylim([0 100])
    xlim([0 nIC])
    set(gca,'fontsize',14,'fontweight','bold')
end
set(gcf,'Position',get(0,'Screensize'))
exportgraphics(gcf,fullfile(Figures_path,'CVs','Fig2_SI.png'))
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 3: CVs SPATIAL MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      CVs SPATIAL MAPS ===================')

tmp         = load_untouch_nii(fullfile(Results_path,['NEW_MASK_ZSCORE_ASSIGN_WM_CSF_MASKED.nii.gz']));
grp_mask    = double(tmp.img);
[rr,cc,zz]  = size(grp_mask);

for k = 1:K
    disp(['State' num2str(k)])
    cv_act_mask         = nan*zeros(size(grp_mask));    
    dist_cv_act_mask    = nan*zeros(size(grp_mask));
    for i = 1:rr
        for j = 1:cc
            for l = 1:zz
                idx_IC = grp_mask(i,j,l);
                if idx_IC==0
                    continue
                else
                    cv_act_mask(i,j,l)      =  new_cv2(k,find(slctIC_origIC==idx_IC));
                    dist_cv_act_mask(i,j,l) =  dist_cv_mean(k,find(slctIC_origIC==idx_IC));

                end
            end
        end
    end
    
%     % 3D NIFTI FILE WITH ALL CVs OF 46 ICs TOGETHER
%     NII = create_4D_nii(fullfile(Results_path,'NEW_MASK_ZSCORE_ASSIGN_WM_CSF_MASKED.nii.gz'),cv_act_mask);
%     save_untouch_nii(NII,fullfile(Figures_path,'CVs',['CV_3D_K' num2str(k) '.nii.gz']));
%     % 3D NIFTI FILE WITH ALL distance-CVs OF 46 ICs TOGETHER
%     NII = create_4D_nii(fullfile(Results_path,'NEW_MASK_ZSCORE_ASSIGN_WM_CSF_MASKED.nii.gz'),dist_cv_act_mask);
%     save_untouch_nii(NII,fullfile(Figures_path,'CVs',['DIST_CV_3D_K' num2str(k) '.nii.gz']));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHRONNECTOME METRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = cell(N,1);
for j=1:N
    T{j} = [Vol];
end

Gamma_reshape = (reshape(Gamma,[Vol N K]));
GammaMeanSubj = squeeze(mean(Gamma_reshape,1));
vpath_reshape = (reshape(vpath,[Vol N]))';

disp(' ')
disp('=============      FRACTIONAL OCCUPANCY, LIFETIME, SWITCHING RATE ===================')

FO = getFractionalOccupancy(Gamma,T,options,2);
for k = 1:K
    [h_k(k),p_k(k)] = kstest2(FO(1:44,k),FO(45:end,k),'Alpha',0.05);
end
[h_k_fdr,p_k_fdr] = fdr_bh(p_k,0.05,'dep','yes');

threshold = round(Vol/100);
for ss =1:N
    switchingRate(ss,1) = getSwitchingRate(squeeze(Gamma_reshape(:,ss,:)),T{ss,1},options);
    lifetimes{ss,:}     = getStateLifeTimes(vpath_reshape(ss,:)',T{ss,1},options);
    intervals{ss,:}     = getStateIntervalTimes(vpath_reshape(ss,:)',T{ss,1},options);
    maxFO(ss,1)         = getMaxFractionalOccupancy(squeeze(Gamma_reshape(:,ss,:)),T{ss,1},options);
    
end

for ss =1:N
    for k = 1:length(lifetimes{ss,1})
        mean_lifetimes{ss,k} = mean(lifetimes{ss,1}{1,k},'omitnan');
    end
end
for k = 1:size(mean_lifetimes,2)
    mean_lifetimes_yo(1,k) = mean([mean_lifetimes{1:44,k}],'omitnan');
    mean_lifetimes_yo(2,k) = mean([mean_lifetimes{45:end,k}],'omitnan');
end
for k = 1:K
    [h_k_LT(k),p_k_LT(k)] = kstest2([mean_lifetimes{1:44,k}],[mean_lifetimes{45:end,k}],'Alpha',0.05);
end
[h_k_LT_fdr,p_k_LT_fdr] = fdr_bh(p_k_LT,0.05,'dep','yes');
save(fullfile(Results_path,'Chronnectome','Metrics.mat'),'FO','switchingRate','lifetimes','intervals','maxFO','mean_lifetimes','mean_lifetimes_yo');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4: VITERBI PATHS AND TP MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      Figure 4 ===================')

group = {'YOUNG','OLD'};
Masks = {[1:size(Gamma,1)/2] [(size(Gamma,1)/2)+1:size(Gamma,1)]};
[P,Pi] = getMaskedTransProbMats(X,cell2mat(T),hmm,Masks,Gamma,Xi);

clear tmp
for jj = 1:size(P,2)    
    tmp(jj,:,:) = P{1,jj};
    tmp(jj,eye(K)==1) = 0;
    for k = 1:K
        tmp(jj,k,:) = tmp(jj,k,:) / sum(squeeze(tmp(jj,k,:)));
    end
end

seconds      = [70 210 350 480 630 770 910];
grid_seconds = seconds/TR;

figure
ax(1) = subplot(221);
imagesc(vpath_reshape(1:length(name_subj_young_25_sel),:));
xticks(grid_seconds), xticklabels(seconds);
yticklabels(''), set(gca,'TickLength',[0,0]);
set(gca,'FontSize',10,'FontWeight','bold');
xlabel('Time (s)','FontSize',14,'FontWeight','bold');
ylabel('\fontsize{15} Subjects','FontSize',14,'FontWeight','bold');
title(['\fontsize{15} VITERBI PATH FOR YOUNG SUBJECTS']);
set(gcf,'Position',get(0,'Screensize'));
color = [0 0 1; 0.643 0 0; 0.343 0 0; 0.1010 0.4450 0.933; 0.85 0.20 0.098; 0.5 0.8 0.3];
cbh = colorbar('XTickLabel',{'S1','S2','S3','S4','S5','S6'},'XTick',[1:K],'Ticks',[1.5 2.25 3.1 3.9 4.75 5.5],...
    'TickLength',0,'FontSize',14,'FontWeight','bold');
colormap(ax(1),color);

ax(2) = subplot(222);
disp_state_matrix(squeeze(tmp(1,:,:)),1:K,[min(tmp(:)) max(tmp(:))]);
xticklabels({'S1','S2','S3','S4','S5','S6'}), yticklabels({'S1','S2','S3','S4','S5','S6'});
title('\fontsize{15} TRANSITION PROBABILITY MATRIX FOR YOUNG SUBJECTS');
colormap(ax(2),jet);

ax(3) = subplot(223);
imagesc(vpath_reshape(length(name_subj_young_25_sel)+1:end,:));
xticks(grid_seconds), xticklabels(seconds);
yticklabels(''), set(gca,'TickLength',[0,0]);
set(gca,'FontSize',10,'FontWeight','bold');
xlabel('Time (s)','FontSize',14,'FontWeight','bold');
ylabel('\fontsize{15} Subjects','FontSize',14,'FontWeight','bold');
title(['\fontsize{15} VITERBI PATH FOR OLD SUBJECTS']);
set(gcf,'Position',get(0,'Screensize'));
color = [0 0 1; 0.643 0 0; 0.343 0 0; 0.1010 0.4450 0.933; 0.85 0.20 0.098; 0.5 0.8 0.3];
cbh = colorbar('XTickLabel',{'S1','S2','S3','S4','S5','S6'},'XTick',[1:K],'Ticks',[1.5 2.25 3.1 3.9 4.75 5.5],...
    'TickLength',0,'FontSize',14,'FontWeight','bold');
colormap(ax(3),color)

ax(4) = subplot(224);
disp_state_matrix(squeeze(tmp(2,:,:)),1:K,[min(tmp(:)) max(tmp(:))]);
xticklabels({'S1','S2','S3','S4','S5','S6'}), yticklabels({'S1','S2','S3','S4','S5','S6'});
title('\fontsize{15} TRANSITION PROBABILITY MATRIX FOR OLD SUBJECTS');
colormap(ax(4),jet);
exportgraphics(gcf,fullfile(Figures_path,'TP_and_ViterbiPath','Fig4.png'))
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5 : FRACTIONAL OCCUPANCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      Figure 5 ===================')

figure
disp_individual_results_boxplot_colored(FO,Nyoung,0.18,{'S1','S2','S3','S4','S5','S6'},'State')
title('DISTRIBUTION OF THE STATES FRACTIONAL OCCUPANCY (FO)')
set(gca,'fontsize',16,'fontweight','bold')
set(gcf,'Position',get(0,'Screensize'));
exportgraphics(gcf,fullfile(Figures_path,'Chronnectome_metrics','Fig5.png'))
%close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 6: MEAN LIFETIMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      Figure 6 ===================')

tmp_box_y = [];
tmp_box_o = [];
for k = 1:K
    tmp1 = cell2mat(mean_lifetimes(1:44,k)); tmp1(isnan(tmp1)) = 0;
    tmp2 = cell2mat(mean_lifetimes(45:end,k)); tmp2(isnan(tmp2)) = 0;
    tmp_box_y = [tmp_box_y; tmp1];
    size_box_y(1,k) = size(tmp1,1);
    tmp_box_o = [tmp_box_o; tmp2];
    size_box_o(1,k) = size(tmp2,1);
    clear tmp1 tmp2
end

cum_size_box_y = cumsum(size_box_y);
cum_size_box_o = cumsum(size_box_o);

G = [TR*tmp_box_y(1:cum_size_box_y(1),:)',TR*tmp_box_o(1:cum_size_box_o(1),:)',...
    TR*tmp_box_y(cum_size_box_y(1)+1:cum_size_box_y(2),:)',TR*tmp_box_o(cum_size_box_o(1)+1:cum_size_box_o(2),:)',...
    TR*tmp_box_y(cum_size_box_y(2)+1:cum_size_box_y(3),:)',TR*tmp_box_o(cum_size_box_o(2)+1:cum_size_box_o(3),:)',...
    TR*tmp_box_y(cum_size_box_y(3)+1:cum_size_box_y(4),:)',TR*tmp_box_o(cum_size_box_o(3)+1:cum_size_box_o(4),:)',...
    TR*tmp_box_y(cum_size_box_y(4)+1:cum_size_box_y(5),:)',TR*tmp_box_o(cum_size_box_o(4)+1:cum_size_box_o(5),:)',...
    TR*tmp_box_y(cum_size_box_y(5)+1:cum_size_box_y(6),:)',TR*tmp_box_o(cum_size_box_o(5)+1:cum_size_box_o(6),:)'];

grp = [ones(1,size_box_y(1,1)),2*ones(1,size_box_o(1,1)),...
    3*ones(1,size_box_y(1,2)),4*ones(1,size_box_o(1,2)),...
    5*ones(1,size_box_y(1,3)),6*ones(1,size_box_o(1,3)),...
    7*ones(1,size_box_y(1,4)),8*ones(1,size_box_o(1,4)),...
    9*ones(1,size_box_y(1,5)),10*ones(1,size_box_o(1,5)),...
    11*ones(1,size_box_y(1,6)),12*ones(1,size_box_o(1,6))];


figure
temp    = [1:K];
temp2   = repmat(temp',1,2)+[0 0.2]; 
positions = reshape(temp2',2*K,1);
boxplot(G,grp,'positions',positions,'Width',0.18,'OutlierSize',7,'Symbol','k+')
color = [1 0 0; 0 0 1]; %% first color is the last to be assigned to the boxplot
color = repmat(color,K,1);
h = findobj(gca,'Tag','Box');
for j = 1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5);
end
c = get(gca, 'Children');
hleg1 = legend(c(1:2),'YOUNG','OLD','Orientation','horizontal','NumColumns',3,'Location','North');
legend boxoff, set(gca,'TickLength',[0,0]), set(gcf,'Position',get(0,'Screensize'))
set(gca,'FontSize',15,'FontWeight','bold')
xticks([1.1:1:6.1]), xticklabels({'S1','S2','S3','S4','S5','S6'}), xlabel('State'), ylabel('Time(s)')
title('DISTRIBUTION OF THE STATES LIME TIME (LT)')
set(gca,'fontsize',16,'fontweight','bold')
set(gcf,'Position',get(0,'Screensize'));
ylim([-5 80])
exportgraphics(gcf,fullfile(Figures_path,'Chronnectome_metrics','Fig6.png'))
%close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 7: MEAN ABSOLUTE BOLD ACTIVATION IN FUNCTIONAL DOMAINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      FIGURE 7 ===================')

for k = 1:K
    disp(['State' num2str(k)])
    pos  = 1;
    for cc = 1:size(comp_network_names,1)
        ll = length(comp_network_names{cc,2});
        avg_bold_dom(cc,k) = 100*mean(abs(hmm.state(k).W.Mu_W(pos:pos+ll-1)));
        pos = pos + ll;
    end
end
max_avg_bold_dom = round(max(avg_bold_dom(:)));
min_avg_bold_dom = round(min(avg_bold_dom(:)));

figure
for k = 1:K
    subplot(3,2,k)
    hold on
    imagesc(avg_bold_dom(:,k))
    colormap('pink')
    colorbar
    caxis([min_avg_bold_dom max_avg_bold_dom])
    axis([0.5 1 0.5 size(comp_network_names,1)+0.5])
    yticks(1:size(comp_network_names,1))
    yticklabels(comp_network_names)
    set(gca,'TickLength',[0,0])
    xticklabels('')
    for cc = 1:size(comp_network_names,1)
        plot([0 1],[cc-0.5 cc-0.5],'Color',[0 0 0],'LineWidth',0.5)
    end
    hold off
    set(gca,'fontsize',12,'fontweight','bold')
    set(gcf,'Position',get(0,'Screensize'));
    title(['S' num2str(k)])
end
exportgraphics(gcf,fullfile(Figures_path,'Mu','Fig7.png'))
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 8: FUNCTIONAL CONNECTIVITY (FC) MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      FIGURE 8 ===================')

figure
for k = 1:K
    subplot(2,3,k)
    disp_basic_matrix(hmm.state(k).W.corr,lst_rsn_name,[-1 1])
    title(['S' num2str(k)],'FontSize',16,'FontWeight','bold')
    set(gcf,'Position',get(0,'Screensize'));
    colormap('jet')
    hold on
    pos_x = 0;
    pos_y = 0;    
    for nn = 1 : length(comp_network_names)
        pos_x = pos_x + length(comp_network_names{nn,2});
        pos_y = pos_y + length(comp_network_names{nn,2});
        if nn < length(comp_network_names)
            plot([pos_x+0.5  pos_x+0.5],[0 nIC+2],'Color',[0 0 0],'LineWidth',0.5)
            plot([0 nIC+2],[pos_y+0.5  pos_y+0.5],'Color',[0 0 0],'LineWidth',0.5)
        end
    end
    hold off
    c = colorbar('Position',...
    [0.925478645066274 0.284061696658098 0.0163475699558175 0.468380462724934],...
    'Ticks',[-1 -0.5 0 0.5 1],...
    'TickLength',0,...
    'TickLabels',{'-1','-0.5','0','0.5','1'},...
    'FontSize',15,'FontWeight','bold');
    hold off
end
exportgraphics(gcf,fullfile(Figures_path,'FC','Fig8.png'))
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAPH METRICS: strengths, weighted local efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      GRAPH METRICS ===================')

s_thr = 0.2; %Sparsity threshold

for k = 1:K
    disp(['State' num2str(k)])

    W_thr       = threshold_proportional(hmm.state(k).Corr, s_thr);
    bin_corr    = W_thr > 0;
    % Strength
    strength(k,:)     = strengths_und(W_thr);
    % Weighted Local Efficiency
    E_wei_local(k,:)  = efficiency_wei(W_thr,2);

end 

% KS-TEST 
for k = 1:K
    k2 = 1;
    while k2 < K+1
        [h_strength(k,k2), p_strength(k,k2)]    = kstest2(squeeze(strength(k,:))',squeeze(strength(k2,:))','Alpha',0.01);
        [h_E_wei_local(k,k2), p_E_wei_local(k,k2)]   = kstest2(squeeze(E_wei_local(k,:))',squeeze(E_wei_local(k2,:))','Alpha',0.01);
        k2 = k2+1;
    end
end
[h_fdr_strength,p_fdr_strength]         = fdr_bh(p_strength,0.05,'dep','yes');
[h_fdr_E_wei_local,p_fdr_E_wei_local]   = fdr_bh(p_E_wei_local,0.05,'dep','yes');

clear h_strength p_strength  h_E_wei_local p_E_wei_local

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAPH METRICS: MEAN AND STD INTO FUNCTIONAL DOMAINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      GRAPH METRICS INTO FUNCTIONAL DOMAINS ===================')

for k = 1:K
    disp(['State' num2str(k)])
    tmp = 0;
    for c = 1:(size(comp_network_names,1))
        rsn = tmp + length(comp_network_names{c,2});
        % MEAN
        avg_strength_RSN(k,c)     = mean(squeeze(strength(k,tmp+1:rsn)));
        avg_E_wei_local_RSN(k,c)  = mean(squeeze(E_wei_local(k,tmp+1:rsn)));
        % STD
        std_strength_RSN(k,c)    = std(squeeze(strength(k,tmp+1:rsn)));
        std_E_wei_local_RSN(k,c) = std(squeeze(E_wei_local(k,tmp+1:rsn)));
        tmp = rsn;
    end
end

% Within states 
% MEAN
avg_state_strength_RSN(1,:)     = mean(avg_strength_RSN,1);
avg_state_E_wei_local_RSN(1,:)  = mean(avg_E_wei_local_RSN,1);
% STD
std_state_strength_RSN(1,:)    = std(avg_strength_RSN,[],1);
std_state_E_wei_local_RSN(1,:) = std(avg_E_wei_local_RSN,[],1);

% Distance 
for k = 1:K
    dist_strength_RSN(k,:)      = (squeeze(avg_strength_RSN(k,:)) - avg_state_strength_RSN) ./ std_state_strength_RSN;
    dist_E_wei_local_RSN(k,:)   = (squeeze(avg_E_wei_local_RSN(k,:)) - avg_state_E_wei_local_RSN) ./ std_state_E_wei_local_RSN;    
end
dist_E_wei_local_RSN(isnan(dist_E_wei_local_RSN))= 0;

save(fullfile(Results_path,'Graph_metrics','Graph_metrics.mat'),'strength','E_wei_local')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 9: GRAPH METRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      Figure 9 ===================')

figure
subplot(221)
disp_group_boxplot_colored(strength,s_thr)
title('STRENGTH','FontSize',18,'FontWeight','bold')
subplot(223)
disp_group_boxplot_colored(E_wei_local,s_thr)
title('WEIGHTED LOCAL EFFICIENCY','FontSize',18,'FontWeight','bold')
subplot(222)
disp_mean_std_area_colored(dist_strength_RSN,[],s_thr,comp_network_names,4)
title('NORMALIZED STRENGTH','FontSize',18,'FontWeight','bold')
subplot(224)
disp_mean_std_area_colored(dist_E_wei_local_RSN,[],s_thr,comp_network_names,4)
title('NORMALIZED WEIGHTED LOCAL EFFICIENCY','FontSize',18,'FontWeight','bold')
exportgraphics(gcf,fullfile(Figures_path,'Graph_metrics','Fig9.png'))
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOUVAIN MODULARITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      MODULARITY ===================')

N_real = 50000; % number of repetitions

% for nr = 1:N_real
%     for k = 1:K
%         disp(['State' num2str(k)])
%         W_thr        = threshold_proportional(hmm.state(k).Corr, 0.2);
%         bin_corr     = W_thr > 0;
%         disp(['Number of negative weights:' num2str(sum(W_thr<0,'all'))])
%         [M_louvain(k,nr,:),Q_louvain(k,nr,:)] = community_louvain(W_thr);
%         M_louvain_FC(nr,k,:,:) = display_matrix_modularity2(hmm.state(k).Corr,squeeze(M_louvain(k,nr,:)),lst_rsn_name_IC);
%         n_modul(nr,k,:) = max(squeeze(M_louvain_FC(nr,k,:,:)),[],'all');
%         
%     end
% end
% 
% save(fullfile(Results_path,'Modularity',['Modularity_' num2str(N_real) 'rep.mat']),'M_louvain_FC','Q_louvain','n_modul','-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODULARITY: MATRIX OF MODE + FILTER AT 95%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(Results_path,'Modularity',['Modularity_' num2str(N_real) 'rep.mat']))

M_louvain_FC_bin = double(M_louvain_FC > 0);
Nreal_95         = N_real*95/100;

for k = 1:K
    disp(['State' num2str(k)])

    M_louvain_FC_mode(k,:,:)       = squeeze(mode(squeeze(M_louvain_FC(:,k,:,:)),1));
    tmp = squeeze(M_louvain_FC_bin(:,k,:,:));
    M_louvain_FC_bin_sum(k,:,:)    = sum(tmp,1);    
    M_louvain_FC_filter_95(k,:,:)  = squeeze(M_louvain_FC_bin_sum(k,:,:))>Nreal_95;
    M_louvain_filtered(k,:,:)      =  M_louvain_FC_mode(k,:,:).*M_louvain_FC_filter_95(k,:,:);
    clear tmp
end

% Modularity values
mean_Q_louvain  = mean(Q_louvain,2);
range_Q_louvain = [min(Q_louvain,[],2) max(Q_louvain,[],2)];

save(fullfile(Results_path,'Modularity',['Results_modularity_' num2str(N_real) 'rep.mat']),'M_louvain_FC_mode','M_louvain_FC_filter_95','M_louvain_filtered','mean_Q_louvain','range_Q_louvain','-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 10: MODULARITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      Figure 10 ===================')

map_modules = [1 1 1; 1 0 1; 0 0.4 0.2;1 0 0; 0.5 0 1; 0.2 0 0];

figure
for k = 1:K
    subplot(2,3,k)
    disp_basic_matrix(squeeze(M_louvain_filtered(k,:,:)),lst_rsn_name,[0 5])
    title(['S' num2str(k)],'FontSize',16,'FontWeight','bold')
    set(gcf,'Position',get(0,'Screensize'));
    colormap(map_modules)
    
    hold on
    pos_x = 0;
    pos_y = 0;    
    for nn = 1 : length(comp_network_names)
        pos_x = pos_x + length(comp_network_names{nn,2});
        pos_y = pos_y + length(comp_network_names{nn,2});
        if nn < length(comp_network_names)
            plot([pos_x+0.5  pos_x+0.5],[0 nIC+2],'Color',[0 0 0],'LineWidth',0.5)
            plot([0 nIC+2],[pos_y+0.5  pos_y+0.5],'Color',[0 0 0],'LineWidth',0.5)
        end
    end
    hold off
    c = colorbar('Position',...
    [0.925478645066274 0.284061696658098 0.0163475699558175 0.468380462724934],...
    'Ticks',[0.5 1.25 2.1 2.9 3.75 4.5],...
    'TickLength',0,...
    'TickLabels',{'0','1','2','3','4','5'},...
    'FontSize',15,'FontWeight','bold');
    c.Label.String = 'Community ID';
    hold off
end
exportgraphics(gcf,fullfile(Figures_path,'Modularity','Fig10.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DICE COMPUTED INTRA FUNCTIONAL DOMAINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('=============      INTRA-DICE COMPUTATION ===================')

for k = 1:K
    for l = 1:K
        tmp = 0;
        for c = 1:(size(comp_network_names,1))
            rsn = tmp + length(comp_network_names{c,2});
            
            % DICE
            DICE_modul_FC_RSN_inside(k,l,c) = dice(squeeze(M_louvain_filtered(k,tmp+1:rsn,tmp+1:rsn))>0,squeeze(M_louvain_filtered(l,tmp+1:rsn,tmp+1:rsn))>0);
            tmp = rsn;
        end
    end
end
clear tmp

% How many stds each state's intra-dice is far from the average computed within the functional domain?
state_pair = [12 13 14 15 16 23 24 25 26 34 35 36 45 46 56];

for c = 1:(size(comp_network_names,1))
    
    current      = DICE_modul_FC_RSN_inside(:,:,c);
    tri_dice_tmp = triu(current,1);
    % MEAN
    avg_state_dice_inside_RSN(c,:) = mean(tri_dice_tmp(tri_dice_tmp>0));
    % STD
    std_state_dice_inside_RSN(c,:) = std(tri_dice_tmp(tri_dice_tmp>0));
    % Distance from the others DICE computed at the state pair level in the
    % specific functional domain
    dist_dice_inside_RSN(:,:,c)    = (current - avg_state_dice_inside_RSN(c,:)) / std_state_dice_inside_RSN(c,:);
    
    clear current tri_dice_tmp current_dist
end


for c = 1:(size(comp_network_names,1))
    cont = 1;
    for k = 1:K
        for l = k+1:K
            DICE_modul_FC_RSN_inside_vect(c,cont)   = DICE_modul_FC_RSN_inside(k,l,c);
            dist_dice_RSN_inside_vect(c,cont)       = dist_dice_inside_RSN(k,l,c);
            cont = cont+1;
        end
    end
end

dist_dice_RSN_inside_vect(isnan(dist_dice_RSN_inside_vect)) = 0;

filt_no_AUD_CER_BG  = DICE_modul_FC_RSN_inside_vect([1:2, 4:9],:);
thr_percentile      = prctile(filt_no_AUD_CER_BG(:),10); %10th percentile of the DICEs distribution

figure
imagesc(filt_no_AUD_CER_BG<=thr_percentile)
yticks(1:8),yticklabels({comp_network_names{[1:2, 4:9],1}})
set(gca,'fontsize',12,'fontweight','bold'),
xticks(1:length(state_pair))
xticklabels({state_pair})
xlabel('State pair')
set(gcf,'Position',get(0,'Screensize'))
exportgraphics(gcf,fullfile(Figures_path,'Modularity','INTRA_Dice_inf_10th_perc.png'))
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DICE COMPUTED INTER FUNCTIONAL DOMAINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:K
    for l = 1:K
        tmp = 0;
        for c = 1:(size(comp_network_names,1))
            rsn = tmp + length(comp_network_names{c,2});
            init_pos = rsn +1;
            if c==1
                DICE_modul_FC_RSN_outside(k,l,c) = dice_mod(squeeze(M_louvain_filtered(k,tmp+1:rsn,[init_pos:nIC])>0),squeeze(M_louvain_filtered(l,tmp+1:rsn,[init_pos:nIC]))>0);
            else
                DICE_modul_FC_RSN_outside(k,l,c) = dice_mod(squeeze(M_louvain_filtered(k,tmp+1:rsn,[1:tmp,init_pos:nIC])>0),squeeze(M_louvain_filtered(l,tmp+1:rsn,[1:tmp,init_pos:nIC]))>0);
            end
            tmp = rsn;
        end
    end
end
clear tmp

% How many stds each state's inter-dice is far from the average computed within the functional domain?
for c = 1:(size(comp_network_names,1))
    
    current      = DICE_modul_FC_RSN_outside(:,:,c);
    tri_dice_tmp = triu(current,1);
    % MEAN
    avg_state_dice_outside_RSN(c,:) = mean(tri_dice_tmp(tri_dice_tmp>0));
    % STD
    std_state_dice_outside_RSN(c,:) = std(tri_dice_tmp(tri_dice_tmp>0));
    % Distance from the others DICE computed at the state pair level in the
    % specific functional domain
    dist_dice_outside_RSN(:,:,c)    = (current - avg_state_dice_outside_RSN(c,:)) ./ std_state_dice_outside_RSN(c,:);
    
    clear current tri_dice_tmp current_dist
end

for c = 1:(size(comp_network_names,1))
    cont = 1;
    for k = 1:K
        for l = k+1:K
            DICE_modul_FC_RSN_outside_vect(c,cont)  = DICE_modul_FC_RSN_outside(k,l,c);
            dist_dice_RSN_outside_vect(c,cont)      = dist_dice_outside_RSN(k,l,c);
            cont = cont+1;
        end
    end
end
dist_dice_RSN_outside_vect(isnan(dist_dice_RSN_outside_vect)) = 0;
dist_dice_RSN_outside_vect(isinf(dist_dice_RSN_outside_vect)) = 0;

filt_all = DICE_modul_FC_RSN_outside_vect;
thr_percentile = prctile(filt_all(not(isnan(filt_all))),10,'all');

figure
imagesc(filt_all<=thr_percentile)
yticklabels({comp_network_names{:,1}})
xticks(1:length(state_pair))
xticklabels({state_pair})
xlabel('State pair')
set(gca,'fontsize',12,'fontweight','bold')
set(gcf,'Position',get(0,'Screensize'))
exportgraphics(gcf,fullfile(Figures_path,'Modularity','INTER_Dice_inf_10th_perc.png'))
close all