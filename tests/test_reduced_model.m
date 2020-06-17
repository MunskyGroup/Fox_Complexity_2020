
clear all
close all
% This is a script to test the approximate FIM, where the FIM has be
% reparameterized. 

% Load the FIM for the full model, plot ellipses.
load('full_FIM.mat')
load('parameters')
load('hog_params')

free_parameters = {'k12','k23','k32','k34','k43','kr2','kr3','kr4','g'};
CV_FIM = inv(FIM);
figure()
for i = 1:9
    for j=1:i
        subplot(9,9,(i-1)*9+j)
        if i==j
           continue
        else            
            hold on
            error_ellipse(CV_FIM([j i],[j i]),log([parameters.(free_parameters{j}),parameters.(free_parameters{i})]));
       end
       if i==9
           xlabel(free_parameters{j})
       end
       if j==1
           ylabel(free_parameters{i})
       end
    end
end

%% Compute the reduced FIM with a neat trick...
% get the eigenvalues of the FIM. 
[vec,vals] = eig(CV_FIM);
% get the loadings from the biggest eigenvalue, and use it to
% reparameterize
% stuff to try: project onto full eigenvector, project onto multiple
% eigenvecxtors. 
loadings = vec(:,1);
%loadings(abs(loadings)<1e-1)=0;
z = [0 0 -log(parameters.kr3)/((log(parameters.k34)+log(parameters.k32))^2) -log(parameters.kr3)/((log(parameters.k34)+log(parameters.k32))^2) 0 0 1/(log(parameters.k34)+log(parameters.k32)) 0 0 ];
% make a jacobian..
J = zeros(length(free_parameters)-2,length(free_parameters));
J(1,1) = 1;
J(2,2) = 1;
%J(3,:) = [0 0 -log(parameters.k34)/(log(parameters.k32).^2) 1/log(parameters.k32) 0 0 0 0 0];
J(3,:) = vec(:,1)';
J(4,:) = vec(:,2)';
J(5,5) = 1;
J(6,6) = 1;
J(7,9) = 1;


% CV_repar = inv(FIM_repar);
%% FIM_repar = J*CV_repar*J';
CV_repar = J*CV_FIM*J';
repar_names = {'k12','k23','alpha1','alpha2','k43','kr2','g'};
figure()
for i = 1:length(CV_repar)
    for j=1:i
        subplot(7,7,(i-1)*7+j)
        if i==j
           continue
        else            
            hold on
            error_ellipse(CV_repar([j i],[j i]),[0,0]);
       end
       if i==length(CV_repar)
           xlabel(repar_names{j})
       end
       if j==1
           ylabel(repar_names{i})
       end
    end
end


    
%% Generate a reduced model. Find slope based on eigenvectors of "true
% parameters?"
alpha1 = parameters.k34/parameters.k32;
alpha2 = parameters.kr3/parameters.k32;

N = 500;


% in the reduction, kr3 and k34 are no longer free parameters. 
free_parameters = {'k12','k23','k32','k43','kr2','g'};
%free_parameters = {'k12','k23','k32','k43','kr2','kr4','g'};

t0 = 2.1e2;
experiment_signal = @(t) hp.A*(((1-exp(-hp.r1*max([0,t-t0])))* ... 
    exp(-hp.r2*max([0,t-t0])))/(1+(((1-exp(-hp.r1*max([0,t-t0])))*exp(-hp.r2*max([0,t-t0])))/hp.M)))^hp.eta;
experiment.input = experiment_signal;
% now, compute the FIM for this subset of parameters
tvec = 60*[ 0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];
experiment.Nc = ones(1,length(tvec)); 
%parameter_ids = [1,3,4,6,8,10,11];
parameter_ids = [1,3,4,6,8,11];

[FIM,FIMs] = get_FIM(parameters,experiment,parameter_ids,N,tvec,1,free_parameters);

CV_FIM = inv(FIM);
figure()
for i = 1:7
    for j=1:i
        subplot(7,7,(i-1)*7+j)
        if i==j
           continue
        else            
            hold on
            error_ellipse(CV_FIM([j i],[j i]),log([parameters.(free_parameters{j}),parameters.(free_parameters{i})]));
       end
       if i==7
           xlabel(free_parameters{j})
       end
       if j==1
           ylabel(free_parameters{i})
       end
    end
end




