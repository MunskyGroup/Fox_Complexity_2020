%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a script to test a number of different boundary conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load true parameters, MLE estimates, etc. 
load('parameters')
load('hog_params')

% get the FSP for the free parameters
t0 = 2.1e2;

experiment_signal = @(t) hp.A*(((1-exp(-hp.r1*max([0,t-t0])))* ... 
    exp(-hp.r2*max([0,t-t0])))/(1+(((1-exp(-hp.r1*max([0,t-t0])))*exp(-hp.r2*max([0,t-t0])))/hp.M)))^hp.eta;
experiment.input = experiment_signal;
N = 500;
tvec = 60*[ 0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];


boundary_conditions = {'FSP','reflecting','renormalize'}
plot_time = 8;
for i=1:length(boundary_conditions)
    %% Get the sensitivity
    experiment.Nc = 500;
    free_parameters = {'k12','kr3'};
    parameter_ids = [1,9];
    [S,z,P] = get_sensitivity(parameters,experiment,parameter_ids,N,1,tvec,1,free_parameters,boundary_conditions{i});
    S1 = reshape(S(:,1),length(S(:,1))/16,16)
    S2 = reshape(S(:,2),length(S(:,2))/16,16)
    subplot(1,3,1)
    hold on
    plot(P(:,plot_time))
    set(gca,'yscale','log')
    xlabel('# RNA')
    ylabel('P(m|t=15)')
    subplot(1,3,2)
    hold on
    plot(abs(S1(:,plot_time)))
    set(gca,'yscale','log')
    xlabel('# RNA')
    ylabel('S_{k21}')
    subplot(1,3,3)
    hold on
    plot(abs(S2(:,plot_time)))
    set(gca,'yscale','log')
    xlabel('# RNA')
    ylabel('S_{kr3}')
end
subplot(1,3,1)
legend(boundary_conditions)

