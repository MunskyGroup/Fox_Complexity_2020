clear all
close all
%% Instead of searching over experiment space, search of the space of degradation rates
ndeg= 50;
load('all_salt_FIMs_deg.mat')
load('stl1_parameters')
load('all_tons')

dt = tons(2)-tons(1); 

% all_degs = linspace(.5*parameters.g,2*parameters.g,ndeg);
all_degs =linspace(.05*parameters.g,5*parameters.g,ndeg);
cells_tmp = 50*ones(1,16);
var_tmp = [];
for i=1:ndeg
        var_tmp(i) = sum(1./(F_STL1_deg(:,:,i)*cells_tmp'));
end

% multply by the required prefactor. 
var_tmp = var_tmp*dt*(1/(tons(end)-tons(1)));
% interpolate to find approximate variance of the true parameters. 
a = find(all_degs>parameters.g);
b = find(all_degs<parameters.g);

y = (var_tmp(a(1))+var_tmp(b(end)))/2;

hold on
plot(all_degs,sqrt(var_tmp))
scatter(parameters.g,sqrt(y))
plot([parameters.g-5.9e-5,parameters.g+5.9e-5],[y,y])
set(gca,'xscale','log','fontsize',14)
xlabel('mRNA degradation rate')
ylabel('Uncertainty about environment')


