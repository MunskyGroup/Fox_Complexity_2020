% compare MLE estimates. 

% load true parameters, MLE estimates, etc. 
load('mle_pars')
load('parameters')
load('hog_params')

% get the FSP for the free parameters
t0 = 2.1e2;

experiment_signal = @(t) hp.A*(((1-exp(-hp.r1*max([0,t-t0])))* ... 
    exp(-hp.r2*max([0,t-t0])))/(1+(((1-exp(-hp.r1*max([0,t-t0])))*exp(-hp.r2*max([0,t-t0])))/hp.M)))^hp.eta;
experiment.input = experiment_signal;
N = 300;
tvec = 60*[ 0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];

%% Get the sensitivity
experiment.Nc = 500;
free_parameters = {'k12','kr3'};
parameter_ids = [1,9];
%[S,z,P] = get_sensitivity(parameters,experiment,parameter_ids,N,1,tvec,1,free_parameters);
 

%%
% Make the FIM
FIM = get_FIM(parameters,experiment,parameter_ids,N,tvec,1,free_parameters);
% FIM = zeros(2,2);
% for k=1:length(tvec)
%     
%     z_t = z(:,k);
%     S_t = S((k-1)*N+1:(k-1)*N+N,:);
%     for i=1:2
%         for j=1:2 
%             z_star = z';
%             FIM(i,j) = FIM(i,j)+experiment.Nc*sum(z_t.*S_t(:,i).*S_t(:,j));
%         end
%     end
% end
%zt = z';
% for i=1:2
%     for j=1:2
%         FIM(i,j) = experiment.Nc*sum(z(:).*S(:,i).*S(:,j));
%     end
% end

%% Plotting
%CV_mle = cov(10.^mle_pars);
CV_mle = cov(log(10.^mle_pars));

hold on
%scatter(10.^mle_pars(:,1),10.^mle_pars(:,2));
scatter(log(10.^mle_pars(:,1)),log(10.^mle_pars(:,2)),[],[0/255,0/255,0/255],'filled');
alpha(.5)

%scatter(parameters.k12,parameters.kr3,'r');
scatter(log(parameters.k12),log(parameters.kr3),[],[255/255,207/255,64/255],'filled');

% error_ellipse(inv(FIM),[parameters.k12,parameters.kr3]);
% error_ellipse(CV_mle,[parameters.k12,parameters.kr3]);
error_ellipse(CV_mle,log([parameters.k12,parameters.kr3]));
error_ellipse(inv(FIM),log([parameters.k12,parameters.kr3]));


xlabel('log k_{12}','Fontsize',16)
ylabel('log k_{r3}','Fontsize',16)



