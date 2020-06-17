%% This is a script to get the FIM. 
clear all

load('mat_files/ctt1_parameters')
load('mat_files/hog_params')

% get the FSP for the free parameters
t0 = 2.1e2;
hp.r2 = 0.0038;
experiment_signal = @(t) hp.A*(((1-exp(-hp.r1*max([0,t-t0])))* ... 
    exp(-hp.r2*max([0,t-t0])))/(1+(((1-exp(-hp.r1*max([0,t-t0])))*exp(-hp.r2*max([0,t-t0])))/hp.M)))^hp.eta;
experiment.input = experiment_signal;
N = 300;
free_parameters = {'k12','k23','k32','k34','k43','kr2','kr3','kr4','g'};
parameter_ids = [1,3,4,5,6,8,9,10,11];
tvec = 60*[ 0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];
experiment.Nc = ones(length(tvec),1);
[FIM,FIMs] = get_FIM(parameters,experiment,parameter_ids,N,tvec,1,free_parameters);

CV = inv(FIM);

npars = length(free_parameters);
for i = 1:npars
    for j=1:i
        if i==npars & j == 1
            figure(2)
            hold on
            %scatter(log(10.^mle_pars(:,j)),log(10.^mle_pars(:,i)),[],[0 0 0],'filled','markerfacealpha',.1);
            error_ellipse(CV([j i],[j i]),log([parameters.(free_parameters{j}),parameters.(free_parameters{i})]));
            xlabel(free_parameters{j})
            ylabel(free_parameters{i})
        end 
        figure(1)
        subplot(npars,npars,(i-1)*npars+j)

        if i==j
            continue
        else
            
            hold on
            error_ellipse(CV([j i],[j i]),log([parameters.(free_parameters{j}),parameters.(free_parameters{i})]));   
        end
       if i==npars
          xlabel(free_parameters{j})
       end
       if j==1
           ylabel(free_parameters{i})
       end
    end
end

