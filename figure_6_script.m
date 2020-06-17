clear all 
close all 

%% Verify optimal experiment design for salt sensing. 
addpath('mat_files','solvers','utils')
tons_long = load('all_tons_long.mat');
load('all_tons.mat')
rng(1);
a = tons(1); b = tons(end); 
n_trials = 50;
nc_total = 1000;
% r1 = rand(1,16);
% r1 = r1./sum(r1);
% 
% r2 = rand(1,16);
% r2 = r2./sum(r2);

tvec = 60*[ 0 1 2 4 6 8 10 15 20 25 30 35 40 45 50 55];
tmin_ind = 6;
tmax_ind = 11;
% unif_same = z>0;
% unif_same = unif_same/sum(unif_same);
unif_tmin_max = zeros(1,length(tvec));
unif_tmin_max(tmin_ind:tmax_ind) = floor(1000*(1/length(unif_tmin_max(tmin_ind:tmax_ind)))); 
% designs = floor(1000.*[z;unif_same;unif_tmin_max;1/16*ones(1,16);r1;r2]);
load all_one_gene_designs
designs(3,:) = unif_tmin_max;


stl1_pdf = load('all_fsp_solutions_stl1_data.mat');
ctt1_pdf = load('all_fsp_solutions_ctt1_data.mat');

stl1_test = load('all_fsp_solutions_stl1_new.mat');
ctt1_test = load('all_fsp_solutions_ctt1_new.mat');
designs = designs(1:5,:);
% add the ctt1 design 
load ctt1_only_opt_design
designs(end+1,:) = floor(z*1000);
load opt_design_2G
twoG = floor(z./sum(z(:))*1000);
designs(end+1,:) = sum(twoG);
ton_store = zeros(length(tons),size(designs,1),n_trials);

%% Optimize the salt using a grid search. this is very slow; results automatically loaded in next section, so not 
% necessary to run. 
designs(3,12)=1;
for ii=1:n_trials
for i=1:length(tons)
    % pick a random salt
%    ind = randi(length(tons));
    ton_i = tons(i);
%     pdf_stl1 = all_fsp_solutions(:,:,ind);
%     pdf_ctt1 = ctt1_solns.all_fsp_solutions(:,:,ind);
    
    pdf_stl1 = stl1_pdf.all_fsp_solutions(:,:,i);
    pdf_ctt1 = ctt1_pdf.all_fsp_solutions(:,:,i);
    
    for j=1:size(designs,1)
        % get the design 
        tmp_design = designs(j,:);
        % simulate a random data set with the given design. 
        if j<6
            pdf = pdf_stl1;
            all_fsp_solutions_tmp = stl1_test.all_fsp_solutions;
        elseif j==6
            pdf = pdf_ctt1;
            all_fsp_solutions_tmp = ctt1_test.all_fsp_solutions;
        else
            % make a pdf with some solutions according to STL1 and some
            % according to CTT1 
            pdf = zeros(size(pdf_stl1));
            all_fsp_solutions_tmp = zeros(size(stl1_test.all_fsp_solutions));
            for k=1:size(pdf,1)
                if twoG(1,k)>0
                    pdf(k,:) = pdf_stl1(k,:);
                    all_fsp_solutions_tmp(k,:,:) = stl1_test.all_fsp_solutions(k,:,:); 
                elseif twoG(2,k)>0
                    pdf(k,:) = pdf_ctt1(k,:);
                    all_fsp_solutions_tmp(k,:,:) = ctt1_test.all_fsp_solutions(k,:,:); 
                end
            end
        end
        
        for k=1:size(pdf,1)
            if tmp_design(k)>0
                sum(pdf(k,:));
                tmp_data(k,:) = sample_d( pdf(k,:), tmp_design(k));
            else
                tmp_data(k,:) = zeros(size(pdf(k,:)));
            end
        end
        
        % find the likelihood of each data set.
        tmp_lhoods = zeros(1,size(all_fsp_solutions_tmp,3));
        for m = 1:size(ctt1_test.all_fsp_solutions,3)
            tmp_lhoods(m) = likelihood(tmp_data,all_fsp_solutions_tmp(:,:,m));
        end
        
        % store the maximum likelihood t_off from the sets. 
        [best,ind2] = max(tmp_lhoods);
        ind2
        % ton_store(i,j,:) = [ton_i tons(ind2)];
        ton_store(i,j,ii) = tons_long.tons(ind2);
       
    end
end
end

%% do some plotting
load all_salt_FIMs
dt = tons(2)-tons(1); 
load('ton_store.mat')
designs(3,12)=1;
designs(3,13)=1;

FI = [];
mse = [];

for i=1:size(designs,1)
    k = size(designs,1)-(i-1);
    figure(100);
    hold on
    scatter(ton_store(:,k,1),ton_store(:,k,2),100,'filled')
    alpha(.5)
    display(['MSE for design ' num2str(i)]) 
    display(['correlation for design ' num2str(i)]) 
    coeff(i) = corr(ton_store(:,i,2),ton_store(:,i,1));
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
    
    l = 1000; b = 1000; w = 275; h=100;
    figure('Position', [l b w h])
    bar(designs(i,:),'k')
    xticks([1:16])
    set(gca,'fontsize',10,'xticklabels',tvec/60)
    legend(['design ' num2str(i)])   
    xlabel('time','fontsize',10)
    ylabel('# of cells','fontsize',10)
    %saveas(gcf,['../../figures/design_bar_' num2str(i)],'epsc')
end
q = sum( (tons'- ton_store).^2 , 3)/(50*64);
mse = sum(q,1);
FI = FI*dt*(1/(tons(end)-tons(1)));
FI = FI([1 2 4 5 6 7]);
mse = mse([1 2 4 5 6 7]); 
designs = designs([1 2 4 5 6 7],:);
% sort the designs by the most informative experiment
[sorted_FI,inds] = sort(FI(1:end),'descend');
sorted_designs = designs(inds,:); 
sorted_mse = mse(inds);

% xlabel('t_{on} true')
% ylabel('t_{on} MLE')
% legend('Design 7','Design 6','Design 5','Design 4', 'Design 3', 'Design 2', 'Design 1')

figure();
clf;
bar(sorted_designs,'k')
xticklabels({'Design 1','Design 2', 'Design 3', 'Design 4', 'Design 5','Design 6','Design 7'})

figure()
yyaxis left
plot(sqrt(sorted_FI),'linewidth',2);
ylabel('FI^{-1}')
xlabel('Experiment ID')

yyaxis right
plot(sqrt(sorted_mse),'linewidth',2);
ylabel('MSE')
xticks(1:5)
xlim([0 6])

figure()
hold on
plot(sqrt(sorted_FI),'linewidth',2);
plot(sqrt(sorted_mse),'linewidth',2);
ylabel('Uncertainty')
xlabel('Experiment ID')

xticks(1:7)
xlim([0 8])

%%
figure()
bar(sqrt(flipud([sorted_FI;sorted_mse]')))
view(90,90)
set(gca,'fontsize',14)



%% Make a plot of the FSP solutions. 
figure
times = [9,10,11,12];
colormap cool
cmap = colormap
for i=1:64
    for j=1:4
    subplot(1,4,j)
    hold on
    plot(ctt1_pdf.all_fsp_solutions(times(j),1:150,i),'color',cmap(i,:))
    ylim([0 0.1])
    title(num2str(tvec(times(j))/60))
    end
end
colorbar('yticklabel',num2cell(round(tons/60.0,2)))


%% Plot a gaussian with opt single-gene design and two gene design.
figure(7);
sigma1 = sorted_FI(end);
sigma2 = sorted_FI(end-1);
sigma3 = sorted_FI(end-1);
ndist = @(x,sigma) ( 1/sqrt(2*sigma) )*exp(-.5*x.^2/(sigma));
x = linspace(-700,700,100)
y1 = ndist(x,sigma1);
y2 = ndist(x,sigma2);
y3 = ndist(x,sigma3);

plot(x,y1);
hold on
plot(x,y2)
plot(x,y3)



%% 
figure();
processed = sorted_designs./sum(sorted_designs,2);
%processed(processed ==0) = 1;
figure
colormap('bone'); imagesc(processed,[-.1,.3]); colorbar
caxis([0 .25])


function l = likelihood(data,marginal)
    marginal=marginal(:,1:size(data,2));
    nz_marginal = marginal>0; % positive (nonzero) probabilities. 
    l = sum(sum(data(nz_marginal).*log10(marginal(nz_marginal))));
end

