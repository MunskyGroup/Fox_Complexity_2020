%%% 
%%  Test the full FSP and compare to the time-scaled FSP. 
load('parameters')
load('hog_params')

% add to the parameter vector:
parameters.k22 = parameters.k23*(parameters.k32/(parameters.k32+parameters.k34));
parameters.k24 = parameters.k23*(parameters.k34/(parameters.k32+parameters.k34));
parameters.k44 = parameters.k43*(parameters.k34/(parameters.k32+parameters.k34));
parameters.k42 = parameters.k43*(parameters.k32/(parameters.k32+parameters.k34));
parameters.kb = parameters.kr3/(parameters.k34+parameters.k43);
parameters.p = parameters.kb / (parameters.kb+1); 
%parameters.kb =(parameters.kr3/(parameters.k32+parameters.k34) + 1)^-1;



t0 = 2.1e2;

experiment_signal = @(t) hp.A*(((1-exp(-hp.r1*max([0,t-t0])))* ... 
    exp(-hp.r2*max([0,t-t0])))/(1+(((1-exp(-hp.r1*max([0,t-t0])))*exp(-hp.r2*max([0,t-t0])))/hp.M)))^hp.eta;
experiment.input = experiment_signal;
N = 250;
tvec = 60*[ 0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];

%% get the reduced model
N = 350;
[A_tv,A_constant,k21,A,A4] = get_A_red(parameters,experiment,N);
p0 = zeros(3*N ,1);
p0(1) = 1.0;
D = 3;

%% solve
ODE = @(t,x) A(t)*x;
[toutred,poutred] = ode23s(ODE,tvec,p0);


%% 
figure()
marginal=zeros(16,N);

for i=1:16
%    signal(i) = experiment_signal(tvec(i));
%    k21g(i) = k21(tvec(i));
%     
    marginal(i,:) = sum(reshape(poutred(i,:),D,N));
    subplot(2,8,i)
    %figure()
    plot(marginal(i,:),'k','Linewidth',4)
    xlim([0,100])
    ylim([0,0.04])
    legend(['t=' num2str(tvec(i)/60)])
    
end

%% Get the full model
N=150;
A  = get_A(parameters,experiment,N);
p0 = zeros(4*N ,1);
p0(1) = 1.0;
D=4;
%% solve
ODE = @(t,x) A(t)*x;
[tout,pout] = ode23s(ODE,tvec,p0);


%% 
marginal=zeros(16,N);

for i=1:16
%    signal(i) = experiment_signal(tvec(i));
%    k21g(i) = k21(tvec(i));
%     
    marginal(i,:) = sum(reshape(pout(i,:),D,N));
    subplot(2,8,i)
    hold on

    %figure()
    plot(marginal(i,:),'r--','Linewidth',2)
    xlim([0,100])
    ylim([0,0.04])
    legend(['t=' num2str(tvec(i)/60)])
    
end
% 
% %% Get the full model
% % rescale some parameters
% parameters.k32 = 1000*parameters.k32;
% parameters.k34 = 1000*parameters.k34;
% parameters.kr3 = 1000*parameters.kr3;
% 
% A  = get_A(parameters,experiment,N);
% p0 = zeros(4*N ,1);
% p0(1) = 1.0;
% D=4;
% %% solve
% ODE = @(t,x) A(t)*x;
% [tout,pout] = ode23s(ODE,tvec,p0);
% 
% 
% %% 
% marginal=zeros(16,N);
% 
% for i=1:16
% %    signal(i) = experiment_signal(tvec(i));
% %    k21g(i) = k21(tvec(i));
% %     
%     marginal(i,:) = sum(reshape(pout(i,:),D,N));
%     subplot(2,8,i)
%     hold on
% 
%     %figure()
%     plot(marginal(i,:),'m:','Linewidth',2)
%     xlim([0,100])
%     ylim([0,0.04])
%     legend(['t=' num2str(tvec(i)/60)])
%     
% end

% figure()
% subplot(1,2,1)
% plot(tvec,signal)
% subplot(1,2,2)
% plot(tvec,k21g)

%% Try getting the sensitivity. 
% free_parameters = [5,9];
% experiment.Nc = 1000;
% [S,z,P] = get_sensitivity(parameters,experiment,free_parameters,N,1,tvec);
% %%
% figure()
% for i=1:16
%     subplot(4,4,i)
%     plot(P(:,i))
%     xlim([0,100])
%     ylim([0,0.04])
%     legend(['t=' num2str(tvec(i)/60)])
%     
% end
% 
% %%
% % Make the FIM
% for i=1:2
%     for j=1:2
%         pp = P';
%         FIM(i,j) = sum(pp(:).*S(:,i).*S(:,j));
%     end
% end
% Plot it



