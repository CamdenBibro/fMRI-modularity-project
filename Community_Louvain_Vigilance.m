close all
clear all

addpath Z:\shared_toolboxes\Derek_functions\
addpath Z:\shared_toolboxes\MatlabProgressBar\
addpath Z:\shared_toolboxes\2019_03_03_BCT\
load("Z:\Camden\fMRI_data\Vigilance_max_min\con_data_max.mat");
load("Z:\Camden\fMRI_data\Vigilance_max_min\pat_data_max.mat");
max_vig_data = permute(con_data_max, [2,3,1]);
pats_2 = permute(pat_data_max, [2,3,1]);
load("Z:\Camden\fMRI_data\Vigilance_max_min\con_data_low.mat");
load("Z:\Camden\fMRI_data\Vigilance_max_min\pat_data_low.mat");
min_vig_data = permute(con_data, [2,3,1]);
pats_1 = permute(pat_data, [2,3,1]);
 
allPatsConsMin = nan(46+51,129,129); 
allPatsConsMax = nan(97,129,129);
allPatsConsMin(1:46,:,:) = con_data; 
allPatsConsMin(47:97,:,:) = pat_data;
allPatsConsMax(1:46,:,:) = con_data_max; 
allPatsConsMax(47:97,:,:) = pat_data_max;
k = isnan(allPatsConsMax(:,1,2));

load("Z:\000_Data\fMRI\spm8_new_preprocessed_data\data\parcellated_and_connectivity_fmri_data\AAN_Bfore_BNST_CONN_2minWindow_2sStride\pat20\pat20_parcellated.mat")
network_names = fMRI_struct.region_names_ic;
net_idx = fMRI_struct.node_network_ids;
[~,net_sorted_idx] = sort(net_idx);
x = 2; 

n = size(max_vig_data,1);
home_bases = ["HG_contra" "HG_ipsi" "aSMG_ipsi" "aSMG_contra" "PP_ipsi" "PP_contra" "CO_ipsi" "CO_contra" "IC_ipsi" "IC_contra" "SMA_ipsi" "SMA_contra"];
targets = ["IFGtri_ipsi" "IFGoper_ipsi" "MedFC_ipsi" "FOrb_ipsi" "toMTG_ipsi" "pSMG_ipsi"];

IDX_home_base = zeros(1,12);
IDX_target = zeros(1,6);
for i = 1:12
    IDX_home_base(:,i) = find(contains(network_names,home_bases(i)));
end
for ii = 1:6
    IDX_target(:,ii) = find(contains(network_names, targets(ii)));
end
%%
% STORING PERMUTATION VALUES
is_target_at_home_question_mark = zeros(1000,1);
M3_original_order = zeros(1000,n); 
M3_original_order_tmp = zeros(1000,n); 


%% MAXIMUM VIGILANCE :: PURMUTATIONS TESTING

iter = 15; % start with 10 as test

n = 129;
pat_data = max_vig_data; % make smaller dataset for testing
[pat] = size(pat_data,3); 

K2 = zeros(n,pat);
gamma1 = 1.1;            % Set GAMMA value

% louvain %
for p = progress(1:pat)
    X = squeeze(pat_data(:,:,p));
    X(isnan(X)) = 0;
    X(X<0) = 0;
    X_scaled = weight_conversion(X, 'normalize');
%       rows_purmute = randperm(n);
   % columns_purmute = randperm(n);
%       X_perm = X_scaled(rows_purmute,rows_purmute); % TRYING ONLY ROW PERMUTATION: Row and column blows up # of communities detected. 
    X_scaled(1:n+1:end) = 0;                        % preserve diagonal
    [M1, Q1a] = consensus_community_louvain_with_finetuning(X_scaled, gamma1); 
%    [~, original_row_idx] = sort(rows_purmute);
%     M1_original = M1(original_row_idx);
    K2(:, p) = M1; 
end
%% RUN PERMUTATIONS

permed_allpat = zeros(1000,n,size(pat_data,3));
for pat = 1:size(pat_data,3)
    for i = 1:1000
        rows_permute = randperm(n);
        permed_allpat(i,:,pat) = K2(rows_permute,pat);
    end
end

k = permed_allpat(:,IDX_home_base,:);
q = permed_allpat(:,IDX_target,:);


%% 
k1 = squeeze(q(:,:,1));


tmp = ones(1000,12,46);
tmp(:,:,1) = tmp(:,:,1)*2;

k1 = permute(q,[1,3,2]);
%k1(:,1,:) = k1(:,1,:)*2;
k2 = reshape(k1,[46000,12]);





%     % association mat across patients %
%     W2 = zeros(size(X));
%     for i = 1:p
%         KK2 = ((squeeze(K2(:,i))) == (squeeze(K2(:,i)))');
%         W2 = W2 + KK2; 
%     end
%     W2 = W2 / p;
%     W3 = weight_conversion(W2, 'normalize');
%      
%     % louvain on association matrix %
%     rows_purmute = randperm(n);
%     [~, original_row_idx] = sort(rows_purmute);
%     W3_perm = W3(rows_purmute, rows_purmute);
%     X_perm(1:n+1:end) = 0;  % preserve diagonal
% 
%   % [M2_max, q1_pats] = consensus_community_louvain_with_finetuning(W2, gamma1);
%     [M3, ~] = consensus_community_louvain_with_finetuning(W3_perm, gamma1); 
%     
%     M3_original_order_tmp(pp,:) = M3(original_row_idx); 
    
% 
% for i = 1:10
%     M3_original_order_tmp(i,IDX_target)
%     M3_original_order_tmp(i,IDX_home_base)
% end

%% 
[~, NMI] = partition_distance(K2');
NMI(NMI==1) = 0;
%figure, imagesc(NMI)
figure, plot(mean(NMI)); title("MIN VIG: mean NMI")

M2_con_max = M2_max == M2_max';
M3_con = M3 == M3';

reordered_NMI = zeros(size(NMI));
reordered_M2_max = M2_max == M2_max'; 
reordered_W3 = zeros(size(NMI)); 
reordered_W2 = zeros(size(NMI));
reordered_M3 = zeros(size(NMI));
[~,idx] = sort(M2_max);
for i = 1:129
    for j = 1:129
        reordered_M2_max(j,i) = M2_con_max(idx(j), idx(i)); 
        reordered_NMI(j,i) = NMI(idx(j), idx(i)); 
        reordered_W3(j,i) = W3(idx(j), idx(i));
        reordered_W2(j,i) = W2(idx(j), idx(i));
        reordered_M3(j,i) = M3_con(idx(j), idx(i));
    end
end
% for i = 1:129
%     for j = 1:129
%         reordered_M2(j,i) = M2_con(net_sorted_idx(j), net_sorted_idx(i)); 
%         reordered_NMI(j,i) = NMI(net_sorted_idx(j), net_sorted_idx(i)); 
%         reordered_W3(j,i) = W3(net_sorted_idx(j), net_sorted_idx(i));
%         reordered_W2(j,i) = W2(net_sorted_idx(j), net_sorted_idx(i));
%         reordered_M3(j,i) = M3_con(net_sorted_idx(j), net_sorted_idx(i));
%     end
% end
figure, imagesc(reordered_M2_max); title("MAX VIG: reordered mod outputs (not normalized)"); xticks([0:129]); xticklabels(network_names(idx)); yticks([0:129]); yticklabels(network_names(idx))
ax = gca;
ax.FontSize = 8;
%figure, imagesc(reordered_M3); title("reordered mod outputs (normalized)")
figure, imagesc(reordered_NMI, "CDataMapping","scaled"); title("max VIG: reordered NMI across nodes"); colorbar
figure, imagesc(reordered_W2, "CDataMapping","scaled"); title("max VIG: reordered association mat from K3 (not normalized)"); colorbar
%figure, imagesc(reordered_W3); title("reordered association mat from K3 (normalized)")

%% Control Louvain %%
n = 129;
pat_data = min_vig_data; % make smaller dataset for testing
[pat] = size(pat_data,3); 

K2 = zeros(n,pat);
gamma1 = 1.1;

% louvain %
for p = progress(1:pat)
    X = squeeze(pat_data(:,:,p));
    X(isnan(X)) = 0;
    X(X<0) = 0;
    X_scaled = weight_conversion(X, 'normalize');
    [M1, Q1a] = consensus_community_louvain_with_finetuning(X_scaled, gamma1); 
    K2(:, p) = M1; 
end

% association mat across patients %
W2 = zeros(size(X));
for i = 1:p
    KK2 = ((squeeze(K2(:,i))) == (squeeze(K2(:,i)))');
    W2 = W2 + KK2; 
end
 W2 = W2 / p;
 W3 = weight_conversion(W2, 'normalize');

 % louvain on association matrix %
[M2_min, q1_con] = consensus_community_louvain_with_finetuning(W2, gamma1);
[M3, ~] = consensus_community_louvain_with_finetuning(W3, gamma1); 

%% 
[~, NMI] = partition_distance(K2');
NMI(NMI==1) = 0;
figure, imagesc(NMI)
figure, plot(mean(NMI)); title("controls mean NMI")

M2_con_min = M2_min == M2_min';
M3_con = M3 == M3';

reordered_NMI = zeros(size(NMI));
reordered_M2_min = M2_min == M2_min'; 
reordered_W3 = zeros(size(NMI)); 
reordered_W2 = zeros(size(NMI));
reordered_M3 = zeros(size(NMI));
[~,idx] = sort(M2_min);
for i = 1:129
    for j = 1:129
        reordered_M2_min(j,i) = M2_con_min(idx(j), idx(i)); 
        reordered_NMI(j,i) = NMI(idx(j), idx(i)); 
        reordered_W3(j,i) = W3(idx(j), idx(i));
        reordered_W2(j,i) = W2(idx(j), idx(i));
        reordered_M3(j,i) = M3_con(idx(j), idx(i));
    end
end
% for i = 1:129
%     for j = 1:129
%         reordered_M2(j,i) = M2_con(net_sorted_idx(j), net_sorted_idx(i)); 
%         reordered_NMI(j,i) = NMI(net_sorted_idx(j), net_sorted_idx(i)); 
%         reordered_W3(j,i) = W3(net_sorted_idx(j), net_sorted_idx(i));
%         reordered_W2(j,i) = W2(net_sorted_idx(j), net_sorted_idx(i));
%         reordered_M3(j,i) = M3_con(net_sorted_idx(j), net_sorted_idx(i));
%     end
% end

figure, imagesc(reordered_M2_min); title("MIN VIG: reordered mod outputs (not normalized)"); xticks([0:129]); xticklabels(network_names(idx)); yticks([0:129]); yticklabels(network_names(idx))
ax = gca;
ax.FontSize = 8;%figure, imagesc(reordered_M3); title("reordered mod outputs (normalized)")
figure, imagesc(reordered_NMI, "CDataMapping","scaled");colorbar; title("MIN VIG: reordered NMI across nodes")
figure, imagesc(reordered_W2, "CDataMapping","scaled"); colorbar; title("MIN VIG: reordered association mat from K3 (not normalized)")
%figure, imagesc(reordered_W3); title("reordered association mat from K3 (normalized)")
