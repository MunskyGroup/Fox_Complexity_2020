%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a script to generate Figure 1 in the manuscript. 
% For the 0.4M and the 0.2M systems.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modify the path
addpath('mat_files','solvers','utils')

%% Write the HOG-MAPK nuclear localization equations
load('mat_files/hog_params')
conditions = [0.2,0.4]; 
r2_params = [0.0069,0.0038]; 
t0 = 2.1e2;
times = linspace(0,60*55,500); 
figure(1)
hold on 
for i=1:length(conditions)
    hp.r2 = r2_params(i);
    experiment_signal = @(t) (t>t0).*hp.A.*(((1-exp(-hp.r1.*(t-t0))).* ... 
        exp(-hp.r2*(t-t0)))./(1+(((1-exp(-hp.r1.*(t-t0))).*exp(-hp.r2.*(t-t0)))./hp.M))).^hp.eta;
    signal = experiment_signal(times);
    plot(times/60,signal)
end
xlabel('time [min]')
ylabel('Nuclear kinase')

%% make some plots of the FSP solutions.
% Choose STL1 or CTT1
analysis_type = 'stl1';
switch analysis_type
case 'stl1'
    load('mat_files/stl1_parameters')
case 'ctt1'
    load('mat_files/ctt1_parameters')
end

model.tvec = 60*[ 0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];
model.parameters = parameters;
model.experiment.input = experiment_signal;
model.N = 200; % FSP dimension
[p_rna,p_full] = solve_fsp_tv(model);

%% plot FSP solutions . 
figure(2)
for j=1:length(model.tvec)
    subplot(2,8,j)
    hold on
    plot(p_rna(j,1:150),'color','k','linewidth',2)
    title(['t=' num2str(model.tvec(j)/60)])
    xlim([0 150])
    ylim([0 0.1])
end
subplot(2,8,9)
xlabel('# mRNA')
ylabel('probability')
