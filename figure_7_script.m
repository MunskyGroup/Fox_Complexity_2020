clear all 
close all 

%% Verify optimal experiment design for salt sensing. 

rng(3)
load all_one_gene_designs
load all_fsp_solutions_experimental
n_trials = 10;
nc_total = 100;
% all_salts = linspace(0.1,0.6,2000);
all_salts = linspace(0.01,0.8,2000); 
tvec = 60*[ 0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];
% Ground truth salt
salts = [0.2, 0.4];
replicate = 1; 

stl1_test = load('all_fsp_solutions_experimental.mat');
ctt1_test = load('all_fsp_solutions_experimental_ctt1.mat');

load all_one_gene_designs
designs = designs(1:5,:);
% add the ctt1 design 
load ctt1_only_opt_design
designs(end+1,:) = floor(z*1000);

load opt_design_2G
twoG = floor(z./sum(z(:))*1000);
designs(end+1,:) = sum(twoG);

%%
%for i=1:n_trials
designs(3,12)=1;
pdf_stl1 = stl1_test.all_fsp_solutions;
pdf_ctt1 = ctt1_test.all_fsp_solutions;
for kk=1:2
    salt = salts(kk); 
    tic;
    for ii=1:n_trials
        % sample a data set 
        for j=1:size(designs,1)
            j
            tic;
            salt_i = salt;
            
            % get the design 
            tmp_design = designs(j,:);
            experiment.Nc = tmp_design;
            experiment.tvec = tvec;
            % simulate a random data set with the given design. 
            if j<6
                gene_id = 1;
                all_fsp_solutions_tmp(:,:,:) = stl1_test.all_fsp_solutions(:,:,:); 
            elseif j==6
                gene_id = 2;
                all_fsp_solutions_tmp(:,:,:) = ctt1_test.all_fsp_solutions(:,:,:);
            else
                % make a pdf with some solutions according to STL1 and some
                % according to CTT1
                pdf = zeros(size(pdf_stl1));
                all_fsp_solutions_tmp = zeros(size(stl1_test.all_fsp_solutions));
                for k=1:size(pdf,1)
                    if twoG(1,k)>0
                        gene_id(k) = 1;
                        all_fsp_solutions_tmp(k,:,:) = stl1_test.all_fsp_solutions(k,:,:); 
                    elseif twoG(2,k)>0
                        gene_id(k) = 2;
                        all_fsp_solutions_tmp(k,:,:) = ctt1_test.all_fsp_solutions(k,:,:);
                    else
                        gene_id(k) = 1;
                    end
                end
            end
            tmp_data = sample_experiment(1, experiment, salt_i, 1, gene_id);
            tmp_lhoods = zeros(1,length(all_salts));
            for i=1:length(all_salts)
                [kk ii j i] 
                % Get the FSP solutions. 
                salt_i = salt;
                
                % find the likelihood of each data set.
                tmp_lhoods(i) = likelihood(tmp_data,all_fsp_solutions_tmp(:,:,i));
            end
            [best,ind2] = max(tmp_lhoods);
            salt_store(ii,j,:,kk) = [salt_i all_salts(ind2) best];
            toc
        end
    end
end


%%
load salt_sensing_exp_verification_v2.mat
load all_tons.mat
FI = [];
mse = [];
load all_salt_FIMs
dt = tons(2)-tons(1); 
load('ton_store_new.mat')
designs(3,12)=1;
designs(3,13)=1;


% Get the FI for each design. 
all_densities = {};
all_XI = {};
all_F = {};
for i=1:size(designs,1)
   k = size(designs,1)-(i-1);
   if i<6
       load('all_salt_FIMs.mat')
   elseif i==6
       load('all_salt_FIMs_CTT1.mat')
       F = F_CTT1;
   elseif i==7
       load('all_salt_FIMs.mat')
       F_STL1 = F;
       load('all_salt_FIMs_CTT1.mat')
       F = zeros(size(F_STL1));
       for k=1:size(F,2)
           if twoG(1,k)>0
               F(:,k) = F_STL1(:,k);
           elseif twoG(2,k)>0
               F(:,k) = F_CTT1(:,k);
           end
       end
   end
   FI(i) = sum(1./(F*designs(i,:)'));
   Q = salt_store(:,i,2,:);
   for jj = 1:2
       [F,XI] = ksdensity( Q(:,:,:,jj) );
       all_XI{i,jj} =  XI;
       all_F{i,jj} = F; 
   end    
    var1(i) = var( salt_store(:,i,2,1) );
    var2(i) = var( salt_store(:,i,2,2) ); 
end



FI = FI*dt*(1/(tons(end)-tons(1)));
FI = FI([1 2 4 5 6 7]);
var1 = var1([1 2 4 5 6 7]); 
var2 = var2([1 2 4 5 6 7]); 
designs = designs([1 2 4 5 6 7],:);
all_XI = all_XI([1 2 4 5 6 7],:);
all_F = all_F([1 2 4 5 6 7],:);

% sort the designs by the most informative experiment
[sorted_FI,inds] = sort(FI(1:end),'descend');
sorted_designs = designs(inds,:); 
sorted_var1 =var1(inds);
sorted_var2 =var2(inds);
all_XI = all_XI(inds,:);
all_F = all_F(inds,:);

figure()
bar(sorted_designs,'k');
cmap = cool;
colors = cmap(1:5:end,:);
count = 1; 
for i=1:size(sorted_designs)
    for jj = 1:2
       figure(jj);
       plot(all_XI{i,jj}, all_F{i,jj}, 'color', colors(count,:), 'linewidth', 2)
       hold on
       xlabel('salt concentration')
       ylabel('frequency')
       xlim([0 .4]+(jj-1)*.3)
       count = count + 1; 
    end
end

figure(1)
legend('Design 1', 'Design 2', 'Design 3', 'Design 4', 'Design 5', 'Design 6', 'Design 7')
figure(2)
legend('Design 1', 'Design 2', 'Design 3', 'Design 4', 'Design 5', 'Design 6', 'Design 7')
%%
figure()
x = linspace(1,10,6);
w = 0.3;
yyaxis left
bar(x,sqrt(sorted_FI),w);
title('0.2 M experiment')
ylabel('predicted stdv (sec)')
xlabel('Experiment ID')

yyaxis right
bar(x+.5,sqrt(sorted_var1),w);
ylabel('measured stdv (M)')
%xticks(1:5)
%xlim([0 6])


figure()
yyaxis left
bar(x,sqrt(sorted_FI),w);
ylabel('predicted stdv (sec)')
xlabel('Experiment ID')
ylim([0,35])

yyaxis right
bar(x+.5,sqrt(sorted_var2),w);
title('0.4 M experiment')
ylabel('measured stdv (M)')
% xticks(1:5)
% xlim([0 6])

figure()
hold on 
scatter(sqrt(sorted_FI),sqrt(sorted_var1),75,'r','filled'); 
scatter(sqrt(sorted_FI),sqrt(sorted_var2),75,'k','filled'); 
xlabel('predicted stdv (sec)')
ylabel('measured stdv (M)')
legend('0.2 M experiments', '0.4 M experiments') 



function l = likelihood(data,marginal)
    if size(data,2) >= size(marginal,2)
        data = data(:,1:size(marginal,2));
    elseif   size(data,2) < size(marginal,2)
        marginal = marginal(:,1:size(data,2));
    end
    nz_marginal = marginal>0; % positive (nonzero) probabilities. 
    l = sum(sum(data(nz_marginal).*log10(marginal(nz_marginal))));
end

