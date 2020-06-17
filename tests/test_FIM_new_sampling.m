% compare MLE estimates. 
clear all
close all

% load true parameters, MLE estimates, etc. 
%load('mle_pars_all_pars')

load('parameters')
load('hog_params')
load('all_fsp_solutions_pars_5x')

% get the FSP for the free parameters
t0 = 2.1e2;

experiment_signal = @(t) hp.A*(((1-exp(-hp.r1*max([0,t-t0])))* ... 
    exp(-hp.r2*max([0,t-t0])))/(1+(((1-exp(-hp.r1*max([0,t-t0])))*exp(-hp.r2*max([0,t-t0])))/hp.M)))^hp.eta;
experiment.input = experiment_signal;
N = 300;
tvec = 60*[ 0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];

%% Get the sensitivity
experiment.Nc = 500;
free_parameters = {'k12','k23','k32','k34','k43','kr2','kr3','kr4','g'};
parameter_ids = [1,3,4,5,6,8,9,10,11];
%%[S,z,P] = get_sensitivity(parameters,experiment,parameter_ids,N,1,tvec,1,free_parameters);
% 
%
%%%
%% Make the FIM
try
    load full_FIM
catch   
    tic
    FIM = get_FIM(parameters,experiment,parameter_ids,N,tvec,0,free_parameters);
    toc
end

%% Make a plot of the covariances. 
figure()
CV_FIM = inv(FIM); 
CV_FIM = flipud(CV_FIM);
CV_FIM = [[CV_FIM zeros(size(CV_FIM,1),1)] ; zeros(1,size(CV_FIM,2)+1)];
pcolor(log(abs(CV_FIM)));
xticks(1.5:1:9.5)
xticklabels(free_parameters)

yticks(1.5:1:9.5)
yticklabels(fliplr(free_parameters))
colorbar

xs = 1.5:1:9.5;
signs = sign(CV_FIM);
for i=1:length(xs)
    for j=1:length(xs)
        if signs(i,j)>0
            note = '+';
        else
            note = '-';
        end
        text(xs(i)-.1,xs(j),note,'color',[1,1,1],'FontSize',24)
    end
end

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
% for each column (data set) find the index of the best params.
[ml,ind_mle] = max(lhood_matrix);
% collect all of the best parameters into a matrix, mle pars
mle_pars = pars_to_test(ind_mle,:);

CV_mle = cov(log(mle_pars));
mu =  mean(log(mle_pars));
% CV FIM
CV_FIM = inv(FIM);

hold on
%scatter(10.^mle_pars(:,1),10.^mle_pars(:,2));
count = 1;
for i = 1:9
    for j=1:i
        subplot(9,9,(i-1)*9+j)
        if i==j
            hist(log(mle_pars(:,i)),25)
        else
            
            hold on
            scatter(log(pars_to_test(:,j)),log(pars_to_test(:,i)),[],'b','filled','markerfacealpha',.1);
            scatter(log(mle_pars(:,j)),log(mle_pars(:,i)),[],[0 0 0],'filled','markerfacealpha',.5);
           
            error_ellipse(CV_FIM([j i],[j i]),log([parameters.(free_parameters{j}),parameters.(free_parameters{i})]));
            error_ellipse(CV_mle([j i],[j i]),log([parameters.(free_parameters{j}),parameters.(free_parameters{i})]));
    
        end
       if i==9
       xlabel(free_parameters{j})
       end
       if j==1
           ylabel(free_parameters{i})
       end
    end
end


%scatter(parameters.k12,parameters.kr3,'r');
%scatter(log(parameters.k12),log(parameters.kr3),'r');

% error_ellipse(inv(FIM),[parameters.k12,parameters.kr3]);
%error_ellipse(CV_mle,[parameters.k12,parameters.kr3]);
%error_ellipse(inv(FIM),log([parameters.k12,parameters.kr3]));
% error_ellipse(CV_mle,log([parameters.k12,parameters.kr3]));
% 
% xlabel('k_{12}','Fontsize',16)
% ylabel('k_{r3}','Fontsize',16)



