%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a script to try different experiments and figure out the optimal
% experiment over some different experiment design spaces. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all 

%% Modify the path
addpath('mat_files','solvers','utils')
load('stl1_parameters') 
load('hog_params')

%% Specify the analysis
analysis = '0.2M'
gene = 'stl1'

switch gene 
    case 'stl1'
    
    switch analysis
        case '0.2M'
            load('mat_files/experimental_FIMs_2M_reduced.mat')
            FIMs = FIMs_2M;
        case '0.4M'
            load('mat_files/experimental_FIMs_4M_reduced.mat')
            FIMs = FIMs_4M;
    end
    
    case 'ctt1'
        
    switch analysis
        case '0.2M'
            load('mat_files/ctt1_fims_all_times_all_params_2M.mat')
            keepers = [1,2,5,6,9];
            for i=1:15
                FIM_i = FIMs{i};
                FIMs{i} = FIM_i(keepers,keepers);
            end
          %  FIMs = FIMs_2M;
        case '0.4M'
            load('mat_files/ctt1_fims_all_times_all_params_4M.mat')
            keepers = [1,2,5,6,9];
            for i=1:15
                FIM_i = FIMs{i};
                FIMs{i} = FIM_i(keepers,keepers);
            end
          %  FIMs = FIMs_4M;
    end
end

%% Greedy algorithm
nsteps = 100000;
% Set the initial state for the algoritm
nc_optimal = zeros(nsteps,length(FIMs));
nc_optimal(2,:) = zeros(1,length(FIMs)); 
% start at random
% tmp = rand(1,length(FIMs));
% tmp = floor(1000*tmp/sum(tmp));
% nc_optimal(1,:) = tmp;

for i=2:nsteps
    E_opt_tmp = zeros(1,length(FIMs));
    % add 1 cell to the best possible time point. 
    for k = 1:length(FIMs)
        % increase the kth entry of the optimal cells by increment
        cells_tmp = nc_optimal(i-1,:);
        cells_tmp(k) = nc_optimal(i-1,k)+1;
        % Get the optimality criteria
        E_opt_tmp(k) = ncell_opt(cells_tmp,FIMs,'subset_volume');
    end
    % update the optimal vector after the ith addition. 
    [~,ind] = max(E_opt_tmp);
    nc_optimal(i,:) = nc_optimal(i-1,:);
    nc_optimal(i,ind) = nc_optimal(i,ind)+1;
end
% Define the optimal design
optimal_design  = nc_optimal(end,:)/sum(nc_optimal(end,:));
switch analysis 
    case '0.2M'
        % Get the E-optimality for the optimal experiment: 
        rep1_nc_2M = [1930 1131 750 996 743 1250 542 822 1150 591 821 385 694 1004 451 1122];
        rep2_nc_2M = [2028  1083 1057 966 884 497 712 575 1004 149 188 560 470 344 84 528]; 
        nc_2M = rep1_nc_2M + rep2_nc_2M; 
        nc_2M = nc_2M(2:end);
        intuitive_design = nc_2M/sum(nc_2M)

    case '0.4M'
        rep1_nc_4M = [927  3391 1043 616 633 457 1817 3512 1924 1554 366 241 543 542 439 466];
        rep2_nc_4M = [368 1417 305 178 737 277 409 305 1187 96 333 340 365 173 505 462];
        nc_4M = rep1_nc_4M+rep2_nc_4M;
        nc_4M = nc_4M(2:end);
        intuitive_design = nc_4M/sum(nc_4M);
    end


[info_expt, expt_FIM] = ncell_opt(floor(intuitive_design*1000),FIMs,'subset_volume');
[info_opt,opt_FIM] = ncell_opt(floor(optimal_design*1000),FIMs,'subset_volume');
disp(['Optimal information: ' num2str(info_opt)])
disp(['Experimental information: ' num2str(info_expt)])

%% Make a plot of the optimal FIM with parameter names and the experimentally measured one.
figure
make_uncertainty_plot(rand(10,5),'scatter_plot',0,'ellipse',0,'crlb',{inv(opt_FIM),inv(expt_FIM)})

%% Make a bar plot comparing the intuitive vs. optimal experiment design. 
figure
bar([info_opt,info_expt])
set(gca,'yscale','log')
ylim([1e8 1e12])

%% Make a plot of how the optimization went
figure
tvec = [1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];
time_labels = {};
for i=1:length(tvec)
    time_labels{i}=[num2str(tvec(i)) ' min'];
end
stairs(nc_optimal./sum(nc_optimal,2))
set(gca,'yscale','log')
xlabel('Iteration')
ylabel('Fraction of cells')

%% Make a plot of the information at each time
crbs = {};
model.free_parameters = {'k12','k23','k32','k34','k43','kr2','kr3','kr4','g'};
model.parameter_ids = [1,3,4,5,6,8,9,10,11];
model.parameters = parameters
keepers = [1,2,5,6,9];
model.free_parameters = model.free_parameters(keepers);
model.parameter_ids = model.parameter_ids(keepers);

diag_data = zeros(5,15); 
for i=1:length(FIMs)
  crbs{i}=inv(FIMs{i});
  diag_data(:,i) = diag(FIMs{i}');
end
figure
free_param_vals = get_param_vec(model.free_parameters,model.parameters);
make_uncertainty_plot(rand(10,5),'scatter_plot',0,'ellipse',1,'crlb',crbs,'true_parameters',free_param_vals)
tvec = 60*[ 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];
    
% find the times that were kept in the optimal design. 
tkeeps = tvec(find(floor(1000*optimal_design)));
for i=1:length(tkeeps)
    FIMs{find(tvec/60==tkeeps(i))}
end 

figure
colormap bone
cmap = colormap;
b = bar(diag_data,'facecolor','flat');
for k = 1:size(diag_data,2)
    b(k).CData = k;
end

legend_str = strsplit(num2str(tvec/60));
for i=1:length(legend_str)
    if ismember(tvec(i)/60,tkeeps)
       legend_str{i}=[legend_str{i} ' min^{**}'];
    else
       legend_str{i}=[legend_str{i} ' min'];
    end
end

legend(legend_str)
model.free_parameters{end} = '\gamma';
xticklabels(model.free_parameters)
ylabel('Diagonal of FI') 
set(gca,'Fontsize',14)
