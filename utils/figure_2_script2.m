%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a script to generate Figures 5+6 in the manuscript. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modify the path
addpath('mat_files','solvers','utils')


% Define the range of the salt values
NaCl = [0,0.6];
% load STL1 parameters
load('../mat_files/parameters')
load('../mat_files/hog_params')
% for .4M

% experimental signal
t0 = 2.1e2;
experiment_signal = @(t) hp.A*(((1-exp(-hp.r1*max([0,t-t0])))* ... 
    exp(-hp.r2*max([0,t-t0])))/(1+(((1-exp(-hp.r1*max([0,t-t0])))*exp(-hp.r2*max([0,t-t0])))/hp.M)))^hp.eta;
model.N = 100;
experiment.input = experiment_signal;
[k21,Acell] = get_A_cell(parameters,experiment,model.N);

tvec = linspace(0,3300,100);
for i=1:length(tvec)
    signal(i) = experiment_signal(tvec(i));
    k21_val(i) = k21(tvec(i));
end

subplot(1,2,1)
plot(tvec/60,signal);
hold on
subplot(1,2,2)
plot(tvec/60,k21_val)
hold on

% find the times where k21>0
pos_k21 = tvec(k21_val>0);
pos_k21_t0 = pos_k21(pos_k21>t0);
ton_2M = pos_k21_t0(1);

% for .4M
hp.r2 = 0.0038;
t0 = 2.1e2;
experiment_signal = @(t) hp.A*(((1-exp(-hp.r1*max([0,t-t0])))* ... 
    exp(-hp.r2*max([0,t-t0])))/(1+(((1-exp(-hp.r1*max([0,t-t0])))*exp(-hp.r2*max([0,t-t0])))/hp.M)))^hp.eta;
experiment.input = experiment_signal;
[k21,Acell] = get_A_cell(parameters,experiment,model.N);

tvec = linspace(0,3300,100);
for i=1:length(tvec)
    signal(i) = experiment_signal(tvec(i));
    k21_val(i) = k21(tvec(i));
end
subplot(1,2,1)
plot(tvec./60,signal);
subplot(1,2,2)
plot(tvec./60,k21_val)
hold on

% find the times where k21>0
pos_k21 = tvec(k21_val>0);
pos_k21_t0 = pos_k21(pos_k21>t0);
ton_4M = pos_k21_t0(1);

% interpolate the tons using 2M and 4M 
m = (ton_4M-ton_2M)/0.2;
b = ton_4M - m*0.4;

% Define a bunch of tons, plot interpolated w/ step assumption
salts = linspace(0.01,.6,500); 
tons = m.*salts + b; 
all_signals = zeros(length(salts),length(tvec)); 
colormap cool
cmap = colormap;
figure
amp = max(k21_val);
for i=1:length(salts)
    all_signals(i,:) = amp*(heaviside(tvec) - heaviside(tvec-t0) + heaviside(tvec-tons(i)));
%     stairs(tvec,all_signals(i,:),'color',cmap(i,:),'linewidth',2);
%     hold on
end
% xlim([0,55])



