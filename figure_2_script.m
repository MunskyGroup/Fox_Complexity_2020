%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a script to generate Figure 1 in the manuscript. 
% For the 0.4M and the 0.2M systems.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all 

%% Modify the path
addpath('mat_files','solvers','utils')
load('stl1_parameters') 

%% Get the model object
model.parameters = parameters;
model.free_parameters = {'k12','k23','k32','k34','k43','kr2','kr3','kr4','g'};
model.parameter_ids = [1,3,4,5,6,8,9,10,11];

% down-select to just those parameters which 
keepers = [1,2,5,6,9];   
model.free_parameters = model.free_parameters(keepers);
model.parameter_ids = model.parameter_ids(keepers);
free_param_vals = get_param_vec(model.free_parameters,model.parameters);

%% Load maximum likelihood results and the corresponding FIM. 
load('experimental_FIMs_2M_reduced.mat')
design='intuitive'
switch design
case 'optimal'
    sim = load('all_mle_pars_opt_design.mat');
    load('optimal_design_stl1_2M.mat')
case 'intuitive'
    sim = load('all_mle_pars_intuitive_exp.mat')
    load('intuitive_design_2M_STL1.mat')
end
[~,full_fim] = ncell_opt(floor(Nc*1000),FIMs_2M,'subset_volume'); % Get the FIM for this experiment design
crlb = inv(full_fim);

%% Plot the uncertainty ellipse for the FIM and for the MLE estimates. 
fig1 = figure;
fig1.Renderer='Painters'; % change the renderer for the scatter plot
make_uncertainty_plot(sim.all_mle_pars,'log_trans',1,'parameter_names',model.free_parameters,'crlb',{inv(full_fim)},'kde',1,'true_parameters',free_param_vals);

%% Get the eigenvalues and eigenvectors for Sigma_mle and FIM^{-1}
[fim_vecs,fim_eigs] = eig(crlb);
[mle_vecs,crb_eigs_sim] = eig(cov(log(sim.all_mle_pars)));

% sort the eigenvalues and eigenvectors
[sorted_fim_eigs,eig_inds] = sort(diag(fim_eigs),'descend');
[sorted_crb_eigs_sim,sim_inds] = sort(diag(crb_eigs_sim),'descend');
sorted_fim_vecs = fim_vecs(:,eig_inds);
sorted_sim_vecs = mle_vecs(:,sim_inds);

% make a bar plot of each eigenvector
figure(1)
phi_good= [];
perp = [];
perp_ang = [];

for i=1:length(sorted_fim_eigs)
    for j=1:length(sorted_fim_eigs)
        perp(i,j) = (sorted_fim_vecs(:,i)' * sorted_sim_vecs(:,j)) / (norm(sorted_fim_vecs(:,i))*norm(sorted_sim_vecs(:,j)));
        perp_ang(i,j) = acosd((sorted_fim_vecs(:,i)' * sorted_sim_vecs(:,j)) / (norm(sorted_fim_vecs(:,i))*norm(sorted_sim_vecs(:,j))));
        % Check if the angles are more than 90 degrees. If so, flip one of
        % the vectors.
        if perp_ang(i,j)>90
            perp_ang(i,j) = acosd((sorted_fim_vecs(:,i)' * -1*sorted_sim_vecs(:,j)) / (norm(sorted_fim_vecs(:,i))*norm(-1*sorted_sim_vecs(:,j))));
        end
        

    end
end

% Get the eigenvectors sorted by their most paralellness
[M,I] = min(abs(1-abs(perp')));
sorted_fim_vecs = sorted_fim_vecs(:,I);
sorted_sim_vecs = sorted_sim_vecs(:,I);

% Plot the eigenvectors
figure()
for i=1:5
    subplot(1,5,i)
    bar([sorted_fim_vecs(:,i),sorted_sim_vecs(:,i)])
end

% Plot the eigenvalues
sorted_sim_vals = sorted_crb_eigs_sim(I); 
figure
bar([sorted_fim_eigs sorted_sim_vals])
for i=1:length(perp)
    ypos = max(sorted_fim_eigs(i),sorted_crb_eigs_sim(i));
    angle = acosd((sorted_fim_vecs(:,i)' * sorted_sim_vecs(:,i)) / (norm(sorted_fim_vecs(:,i))*norm(sorted_sim_vecs(:,i))));  
    if angle>90
        angle = angle-180;
    end
    text(i-.2,ypos+.2*ypos,['\phi=' num2str(round(angle,3,'significant'))])
end
ylim([5e-4,.1])
set(gca,'yscale','log')
