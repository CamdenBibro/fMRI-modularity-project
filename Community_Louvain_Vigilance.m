close all
clear all

addpath Z:\shared_toolboxes\Derek_functions\
addpath Z:\shared_toolboxes\MatlabProgressBar\
addpath Z:\shared_toolboxes\2019_03_03_BCT\
load("Z:\Camden\fMRI_data\Vigilance_max_min\con_data_max.mat");
max_vig_data = permute(con_data_max, [2,3,1]);
load("Z:\Camden\fMRI_data\Vigilance_max_min\con_data_low.mat");
min_vig_data = permute(con_data, [2,3,1]);
load("Z:\000_Data\fMRI\spm8_new_preprocessed_data\data\parcellated_and_connectivity_fmri_data\AAN_Bfore_BNST_CONN_2minWindow_2sStride\pat20\pat20_parcellated.mat")
network_names = fMRI_struct.region_names_ic;
net_idx = fMRI_struct.node_network_ids;
[~,net_sorted_idx] = sort(net_idx);
x = 2; 
% 
% 
% 
% patient_data =  permute(patient_data, [2,3,1]);
% patient_data(:,:,10) = [];        % remove problematic subjects
% control_data = permute(control_data, [2,3,1]); 
% control_data(:,:,[37,5,12]) = []; % remove problematic subjects

%% max vig
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
    scaled_X = weight_conversion(X, 'normalize');
    [M1, Q1a] = consensus_community_louvain_with_finetuning(scaled_X, gamma1); 
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
[M2_max, q1_pats] = consensus_community_louvain_with_finetuning(W2, gamma1);
[M3, ~] = consensus_community_louvain_with_finetuning(W3, gamma1); 

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
    scaled_X = weight_conversion(X, 'normalize');
    [M1, Q1a] = consensus_community_louvain_with_finetuning(scaled_X, gamma1); 
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
