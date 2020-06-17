%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a script to try different experiments and figure out the optimal
% experiment over some different experiment design spaces. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all 

%% Modify the path
addpath('mat_files','solvers','utils')


% Set the rng for reproducibilty of random experiments.
rng(0)

% Pick analysis type: one of '2M', '4M', or 'both'
% generate a bunch of experimental settings. 
analysis_type = '2M';

rep1_nc_2M = [1930 1131 750 996 743 1250 542 822 1150 591 821 385 694 1004 451 1122];
rep2_nc_2M = [2028  1083 1057 966 884 497 712 575 1004 149 188 560 470 344 84 528]; 
nc_2M = rep1_nc_2M + rep2_nc_2M; 

rep1_nc_4M = [927  3391 1043 616 633 457 1817 3512 1924 1554 366 241 543 542 439 466];
rep2_nc_4M = [368 1417 305 178 737 277 409 305 1187 96 333 340 365 173 505 462];
nc_4M = rep1_nc_4M+rep2_nc_4M; 

Nc_total = 1000;
nc_measured(1,:) = nc_2M(2:end);
nc_measured(2,:) = nc_4M(2:end);

n_exp = 1000; 
all_times = 60*[ 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];
for i=1:n_exp
    if strcmp(analysis_type,'both')
        % pick either .2 or .4 M
        all_exp{i}.salt = round(rand)*.2+.2;
        if all_exp{i}.salt==0.4
            ind = 2;
        elseif all_exp{i}.salt==0.2
            ind = 1;
        end
    elseif strcmp(analysis_type,'2M')
        all_exp{i}.salt = 0.2;
        ind = 1;
    elseif strcmp(analysis_type,'4M')
        all_exp{i}.salt = 0.4;
        ind = 2; 
    end
    % pick a random number of times
    nt = randi([2,length(all_times)]);
    
    % pick random times
    all_exp{i}.time_inds = sort(randsample(length(all_times),nt));
    all_exp{i}.times = all_times(all_exp{i}.time_inds);
    % split the total number of cells up
%    while min(nc)>floor(Nc_total/(2*nt))
    nc = diff(sort([0 randi([floor(Nc_total/(2*nt)),Nc_total],1,nt-1) Nc_total]));
    while sum(nc>nc_measured(ind,all_exp{i}.time_inds))>0
        disp('trying again')
        %nc = diff(sort([0 randi([floor(Nc_total/(2*nt)),Nc_total],1,nt-1) Nc_total]));
        nc = diff(sort([0 randi([10,Nc_total],1,nt-1) Nc_total]));
    end
    disp('***********************')
    all_exp{i}.Nc = nc;
end

% Run these experiments through FIM.
% STL1 data structure
%all_FIMs_STL1 = zeros(9,9,length(all_exp));
all_FIMs_STL1 = zeros(5,5,length(all_exp));

all_E_opts_STL1 = zeros(1,length(all_exp));
all_D_opts_STL1 = zeros(1,length(all_exp));
all_g_STL1 = zeros(1,length(all_exp));

% CTT1 data structure
%all_FIMs_CTT1 = zeros(9,9,length(all_exp));
all_FIMs_CTT1 = zeros(5,5,length(all_exp));

all_E_opts_CTT1 = zeros(1,length(all_exp));
all_D_opts_CTT1 = zeros(1,length(all_exp));
all_g_CTT1 = zeros(1,length(all_exp));

inds_2M = [];
inds_4M = [];
all_ntimes =[];
keepers = [1,2,5,6,9];


for i=1:length(all_exp)
    % Make sure number of cells looks good: 
    if sum(all_exp{i}.Nc)~=1000
        display('Number of cells is wrong')
    end
    all_ntimes = [all_ntimes length(all_exp{i}.times)];
    % choose which salt
    if all_exp{i}.salt == .2
        stl1_fims = load('experimental_FIMs_2M.mat');
        ctt1_fims = load('ctt1_fims_all_times_all_params_2M.mat');
        inds_2M = [inds_2M i];
        
        
    elseif all_exp{i}.salt == .4
        stl1_fims = load('experimental_FIMs_4M.mat');
        ctt1_fims = load('ctt1_fims_all_times_all_params_4M.mat');
        inds_4M = [inds_4M i];
    end

    % compute the FIM
    for j = 1:length(all_exp{i}.time_inds)
        
        all_exp{i}.time_inds(j);
        tmpfims_j = stl1_fims.FIMs{all_exp{i}.time_inds(j)};
        all_FIMs_STL1(:,:,i) = all_exp{i}.Nc(j)*tmpfims_j(keepers,keepers);
        tmpfims_j = ctt1_fims.FIMs{all_exp{i}.time_inds(j)};
        all_FIMs_CTT1(:,:,i) = all_exp{i}.Nc(j)*tmpfims_j(keepers,keepers);
        
    end
    E_stl1 = eig(all_FIMs_STL1(:,:,i));
    E_ctt1 = eig(all_FIMs_CTT1(:,:,i));

    all_E_opts_STL1(i) = min(E_stl1);
    all_D_opts_STL1(i) = prod(E_stl1);
    
    all_E_opts_CTT1(i) = min(E_ctt1);
    all_D_opts_CTT1(i) = prod(E_ctt1);
    
end


% save the random experiments. 
save(['mat_files/all_random_experiments_' analysis_type '.mat'],'all_exp')
save(['mat_files/all_FIM_metrics_' analysis_type '.mat'],'all_E_opts_STL1','all_D_opts_STL1')


% Load the required data
results_2M = load('mat_files/all_FIM_metrics_2M.mat');
results_4M = load('mat_files/all_FIM_metrics_4M.mat');

% Load the experimental and optimal information for the system.
opt_design_info_2M = load('mat_files/optimal_experiment_info_0.2.mat');
opt_design_info_4M = load('mat_files/optimal_experiment_info_0.4.mat');

exp_design_info_2M = load('mat_files/intuitive_experiment_info_0.2.mat');
exp_design_info_4M = load('mat_files/intuitive_experiment_info_0.4.mat');

% Make stair plot
figure(1);
hold on
semilogy(sort(results_2M.all_D_opts_STL1,'descend'),'k','linewidth',2)
hold on
semilogy(sort(results_4M.all_D_opts_STL1,'descend'),'g','linewidth',2)


% add horizontal lines for some optimal information values.
semilogy([1,length(results_2M.all_D_opts_STL1)],[opt_design_info_2M.info_opt opt_design_info_2M.info_opt],'k--','linewidth',2)
semilogy([1,length(results_2M.all_D_opts_STL1)],[exp_design_info_2M.info_expt exp_design_info_2M.info_expt],'r--','linewidth',2)

semilogy([1,length(results_4M.all_D_opts_STL1)],[opt_design_info_4M.info_opt opt_design_info_4M.info_opt],'b--','linewidth',2)
semilogy([1,length(results_4M.all_D_opts_STL1)],[exp_design_info_4M.info_expt exp_design_info_4M.info_expt],'m--','linewidth',2)

legend('random 2M','random 4M','optimal 2M','experimental 2M','optimal 4M', 'experimental 4M')
set(gca,'fontsize',14)
xlabel('Experiment ID','Fontsize',18)
ylabel('D-optimality','Fontsize',18)
set(gca,'yscale','log')
% make an added inset
figure(2)
hold on
semilogy(sort(results_2M.all_D_opts_STL1,'descend'),'k','linewidth',2)
hold on
semilogy(sort(results_4M.all_D_opts_STL1,'descend'),'g','linewidth',2)
set(gca,'fontsize',14)
xlabel('Experiment ID','Fontsize',14)
ylabel('D-optimality','Fontsize',14)

% add horizontal lines for some optimal information values.
hold on
semilogy([1,length(results_2M.all_D_opts_STL1)],[opt_design_info_2M.info_opt opt_design_info_2M.info_opt],'k--','linewidth',2)
semilogy([1,length(results_2M.all_D_opts_STL1)],[exp_design_info_2M.info_expt exp_design_info_2M.info_expt],'r--','linewidth',2)

semilogy([1,length(results_4M.all_D_opts_STL1)],[opt_design_info_4M.info_opt opt_design_info_4M.info_opt],'b--','linewidth',2)
semilogy([1,length(results_4M.all_D_opts_STL1)],[exp_design_info_4M.info_expt exp_design_info_4M.info_expt],'m--','linewidth',2)
xlim([1,50])

set(gca,'yscale','log')

legend('random 2M','random 4M','optimal 2M','experimental 2M','optimal 4M', 'experimental 4M')

% make a fancy matrix
% use the 2M FIMs:
load('mat_files/experimental_FIMs_2M_reduced.mat');
load('mat_files/experimental_FIMs_4M_reduced.mat');

info_matrix = zeros(4,2);
info_matrix(1,1) = ncell_opt(opt_design_info_2M.nc_opt,FIMs_2M,'subset_volume');
info_matrix(1,2) = ncell_opt(opt_design_info_2M.nc_opt,FIMs_4M,'subset_volume');

info_matrix(2,1) = ncell_opt(opt_design_info_4M.nc_opt,FIMs_2M,'subset_volume');
info_matrix(2,2) = ncell_opt(opt_design_info_4M.nc_opt,FIMs_4M,'subset_volume');

info_matrix(3,1) = ncell_opt(exp_design_info_2M.nc_1000,FIMs_2M,'subset_volume');
info_matrix(3,2) = ncell_opt(exp_design_info_2M.nc_1000,FIMs_4M,'subset_volume');

info_matrix(4,1) = ncell_opt(exp_design_info_4M.nc_1000,FIMs_2M,'subset_volume');
info_matrix(4,2) = ncell_opt(exp_design_info_4M.nc_1000,FIMs_4M,'subset_volume');

figure(3)
colormap bone
imagesc(log10(info_matrix))
colorbar
set(gca,'fontsize',14)


%% use the STL1 and CTT1 FIMs
load('mat_files/experimental_FIMs_2M_reduced.mat');
load('mat_files/experimental_FIMs_4M_reduced.mat');

ctt1_2M = load('mat_files/ctt1_fims_all_times_all_params_2M.mat');
ctt1_4M = load('mat_files/ctt1_fims_all_times_all_params_4M.mat');

keepers = [1,2,5,6,9];
for i=1:15
    FIM_i = ctt1_2M.FIMs{i};
    ctt1_2M.FIMs{i} = FIM_i(keepers,keepers);
    FIM_i = ctt1_4M.FIMs{i};
    ctt1_4M.FIMs{i} = FIM_i(keepers,keepers);
end


info_matrix = zeros(4,4);
info_matrix(1,1) = ncell_opt(opt_design_info_2M.nc_opt,FIMs_2M,'subset_volume');
info_matrix(1,2) = ncell_opt(opt_design_info_2M.nc_opt,FIMs_4M,'subset_volume');

info_matrix(2,1) = ncell_opt(opt_design_info_4M.nc_opt,FIMs_2M,'subset_volume');
info_matrix(2,2) = ncell_opt(opt_design_info_4M.nc_opt,FIMs_4M,'subset_volume');

info_matrix(3,1) = ncell_opt(exp_design_info_2M.nc_1000,FIMs_2M,'subset_volume');
info_matrix(3,2) = ncell_opt(exp_design_info_2M.nc_1000,FIMs_4M,'subset_volume');

info_matrix(4,1) = ncell_opt(exp_design_info_4M.nc_1000,FIMs_2M,'subset_volume');
info_matrix(4,2) = ncell_opt(exp_design_info_4M.nc_1000,FIMs_4M,'subset_volume');

% Apply STL1 design to CTT1 experiments. 
info_matrix(1,3) = ncell_opt(opt_design_info_2M.nc_opt,ctt1_2M.FIMs,'subset_volume');
info_matrix(1,4) = ncell_opt(opt_design_info_2M.nc_opt,ctt1_4M.FIMs,'subset_volume');

info_matrix(2,3) = ncell_opt(opt_design_info_4M.nc_opt,ctt1_2M.FIMs,'subset_volume');
info_matrix(2,4) = ncell_opt(opt_design_info_4M.nc_opt,ctt1_4M.FIMs,'subset_volume');

info_matrix(3,3) = ncell_opt(exp_design_info_2M.nc_1000,ctt1_2M.FIMs,'subset_volume');
info_matrix(3,4) = ncell_opt(exp_design_info_2M.nc_1000,ctt1_4M.FIMs,'subset_volume');

info_matrix(4,3) = ncell_opt(exp_design_info_4M.nc_1000,ctt1_2M.FIMs,'subset_volume');
info_matrix(4,4) = ncell_opt(exp_design_info_4M.nc_1000,ctt1_4M.FIMs,'subset_volume');

% Apply CTT1 design to STL1 experiments 
load('mat_files/ctt1_optimal_design_2M.mat')
info_matrix(5,1) = ncell_opt(floor(1000*optimal_design),FIMs_2M,'subset_volume');
info_matrix(5,2) = ncell_opt(floor(1000*optimal_design),FIMs_4M,'subset_volume');

info_matrix(5,3) = ncell_opt(floor(1000*optimal_design),ctt1_2M.FIMs,'subset_volume');
info_matrix(5,4) = ncell_opt(floor(1000*optimal_design),ctt1_4M.FIMs,'subset_volume');

% Apply CTT1 design to STL1 experiments 
load('mat_files/ctt1_optimal_design_4M.mat')
info_matrix(6,1) = ncell_opt(floor(1000*optimal_design),FIMs_2M,'subset_volume');
info_matrix(6,2) = ncell_opt(floor(1000*optimal_design),FIMs_4M,'subset_volume');

info_matrix(6,3) = ncell_opt(floor(1000*optimal_design),ctt1_2M.FIMs,'subset_volume');
info_matrix(6,4) = ncell_opt(floor(1000*optimal_design),ctt1_4M.FIMs,'subset_volume');

figure(3)
colormap bone
imagesc(log10(info_matrix))
colorbar
set(gca,'fontsize',14)
