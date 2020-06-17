% test sensitivity
% this is a script to test the sensitivity analysis for the hog model. 


%%  Test the FSP
load('parameters')
load('hog_params')

t0 = 2.1e2;

experiment_signal = @(t) hp.A*(((1-exp(-hp.r1*max([0,t-t0])))* ... 
    exp(-hp.r2*max([0,t-t0])))/(1+(((1-exp(-hp.r1*max([0,t-t0])))*exp(-hp.r2*max([0,t-t0])))/hp.M)))^hp.eta;
experiment.input = experiment_signal;
N = 200;
tvec = 60*[ 0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];
A  = get_A(parameters,experiment,N);
p0 = zeros(4*N ,1);
p0(1) = 1.0;
%% solve the unperturbed FSP
ODE = @(t,x) A(t)*x;
JAC = @(t,x) A(t);
options = odeset('jacobian',JAC);
[tout,pout1] = ode23s(ODE,tvec,p0,options);
% convert to marginal
P = zeros(N,length(tvec));
for i=1:N
    P(i,:) = sum(pout1(:,(i-1)*4+1:(i-1)*4+4)');
end
pout1 = P;

%% perturb parameters, solve the new model.
par_pert = parameters;
% delk12 = 1e-5*parameters.k12;
% par_pert.k12 = parameters.k12+delk12;
delkr3 = 1e-4*parameters.kr3;
par_pert.kr3 = parameters.kr3+delkr3; 
A = get_A(par_pert,experiment,N);
ODE = @(t,x) A(t)*x;
JAC = @(t,x) A(t);
options = odeset('jacobian',JAC);
[tout,pout2] = ode23s(ODE,tvec,p0,options);
% convert to marginal
P = zeros(N,length(tvec));
for i=1:N
    P(i,:) = sum(pout2(:,(i-1)*4+1:(i-1)*4+4)');
end
pout2 = P;
dpdk12 = (pout2-pout1)./delkr3;

%%  get the sensitivity
free_parameters = [9];
experiment.Nc = 1; 
[S,z,P] = get_sensitivity(parameters,experiment,free_parameters,N,1,tvec,0,{'kr3'});

%dpdk12=dpdk12';
%% plot the errors
% S = reshape(S',size(dpdk12(:)));
% er1 = sum(abs(dpdk12-S),2);
% er2 = max(abs(dpdk12-S),[],2);
er1 = sum(abs(dpdk12(:)-S));
er2 = max(abs(dpdk12(:)-S));
% subplot(1,2,1)
% plot(tvec,er1);
% xlabel('time (s)','Fontsize',20)
% ylabel('|S_{FSP} - S_{FD}|_{1}','Fontsize',20)
% subplot(1,2,2)
% plot(tvec,er2);
% xlabel('time (s)','Fontsize',20)
% ylabel('|S_{FSP} - S_{FD}|_{\infty}','Fontsize',20)
% all_er = abs(dpdk12-S)';

figure()
[sorted_S,s_inds] = sort(abs(S));
[sorted_F,f_inds] = sort(abs(dpdk12(:)));
subplot(1,2,1)
hold on
semilogy(sort(abs(S)))
semilogy(sort(abs(dpdk12(:))))
set(gca,'yscale','log')
subplot(1,2,2)
hold on 
semilogy(abs(S))
semilogy(abs(dpdk12(:)))
set(gca,'yscale','log')

figure()
hold on
S = reshape(S,200,16);
for i=1:16
scatter(abs(S(:,i)),abs(dpdk12(:,i)))
end

set(gca,'yscale','log','xscale','log')
plot([min(abs(S)) max(abs(S))],[min(abs(S)) max(abs(S))],'k--','linewidth',4)
xlim([1e-35,1e5])
ylim([1e-35,1e5])


% imagesc(all_er(1:200,:))





% free_parameters = [5,9];
% experiment.Nc = 500;
% [S,z,P] = get_sensitivity(parameters,experiment,free_parameters,N,1,tvec);
