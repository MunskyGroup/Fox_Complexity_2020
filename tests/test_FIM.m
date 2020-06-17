%% This is a script to simulate data sets, fit them, and then compare to the FSP-FIM. 

% Load and generate data
load('parameters')
load('hog_params')

t0 = 2.1e2;

experiment_signal = @(t) hp.A*(((1-exp(-hp.r1*max([0,t-t0])))* ... 
    exp(-hp.r2*max([0,t-t0])))/(1+(((1-exp(-hp.r1*max([0,t-t0])))*exp(-hp.r2*max([0,t-t0])))/hp.M)))^hp.eta;
experiment.input = experiment_signal;
N = 250;
tvec = 60*[ 0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];
% Get the A matrix
A  = get_A(parameters,experiment,N);
[A_tv,A_constant,k21]  = get_A_pw(parameters,experiment,N);
p0 = zeros(4*N ,1);
p0(1) = 1.0;

% Solve
ODE = @(t,x) (A_constant+A_tv*k21(t))*x;
JAC = @(t,x) (A_constant+A_tv*k21(t));
options = odeset('jacobian',JAC);
tic;
[tout,pout] = ode23s(ODE,tvec,p0,options);
solve_time = toc
% get the marginal
marginal=zeros(length(tout),N);
for i=1:length(tout)     
    marginal(i,:) = sum(reshape(pout(i,:),4,N));
end

%% fitting starts here
free_params = {'k12','kr3'};
% get the parameter guess
for i=1:length(free_params)
    pguess(i) = parameters.(free_params{i});
end

n_data = 2;
mle_pars = zeros(n_data,length(free_params));
for i=1:n_data
    % sample the data
    data = zeros(size(marginal));
    for j=1:length(tout)
        di = sample_d(marginal(j,:),500);
        data(j,1:length(di)) = di; 
    end
    clf
    hold on
    plot(data(8,:)./500.0);
    plot(marginal(8,:))
    pause
    % fit the model
%    f = @(x) get_likelihood(10.^x,parameters,free_params,data,experiment,tvec,N);
%    mle_pars(i,:) = fminsearch(f,log10(pguess));

end

   save('mle_pars','mle_pars');
