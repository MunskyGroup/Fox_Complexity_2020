function [] = make_uncertainty_plot(data,varargin)
%%%%%%%%%%%%%%%%%%%
% Make a plot of uncertainties for a given covariance
% data is a matrix with n_samples rows and n_parameters columns
% log_trans is a boolean; if true, take the log of the data, otherwise false. 
% scatter_plot is a boolean; if true, add a scatter plot
% VARARGIN:
% 
%%%%%%%%%%%%%%%%%%%

options = struct('log_trans',0,'scatter_plot',1,'contour_plot',0,'ellipse',1,'parameter_names',{{}},'crlb',[],'kde',0,'true_parameters',[],'colors',[]);

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %
   inpName = lower(pair{1}); %# make case insensitive
   if any(strcmp(inpName,optionNames))
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end

% add parameter names
npars = size(data,2);
if length(options.parameter_names)==0
    for i=1:npars
        options.parameter_names{i} = ['\lambda_' num2str(i)];
    end    
end

if options.log_trans
    data = log(data);
    options.true_parameters = log(options.true_parameters);
end


c = parula;

if length(options.colors)==0
    options.colors=c(randi(64),:);
end
options.colors
mu = mean(data);
CV = cov(data);

hist_color =   [253,233,210]/255;
data_ellipse_color = [249,169,75]/255;
crb_ellipse_colors = {[131,139,197]/255 ,  [0,0,0]};

for i = 1:npars
    for j=1:i
        subplot(npars,npars,(i-1)*npars+j)
        if i==j  
            if options.kde
                hold on
                % do some kernel density estimation
                [f,xi] = ksdensity(data(:,i));                                   
                fill([xi fliplr(xi)],[f zeros(1,length(f))],hist_color);
                if length(options.true_parameters)>0
                    plot([options.true_parameters(i) options.true_parameters(i)],[0,max(f)])
                end
            else
                % or a histogram as default. 
                histogram(data(:,i))
                if length(options.true_parameters)>0
                    hold on
                    N = histcounts(data(:,i));
                    plot([options.true_parameters(i) options.true_parameters(i)],[0,max(N)])
                end
            end
       
        else            
            hold on
            if options.ellipse
                error_ellipse(CV([j i],[j i]),[mu(j),mu(i)],'conf',.9,'color',data_ellipse_color);  
                %error_ellipse(CV([j i],[j i]),[mu(j),mu(i)],'conf',.9,'color',c(40,:));  
            end
            if length(options.crlb)>0
                hold on
                for n=1:length(options.crlb)
                    crlb_n = options.crlb{n};
                    if n>length(crb_ellipse_colors)
                    error_ellipse(crlb_n([j i],[j i]),[mu(j),mu(i)],'conf',.9,'color',crb_ellipse_colors{1});
                    else
                    error_ellipse(crlb_n([j i],[j i]),[mu(j),mu(i)],'conf',.9,'color',crb_ellipse_colors{n});
                    end
                end
            end
                
            if options.scatter_plot
                hold on
                scatter(data(:,j),data(:,i),10,'k','filled')
                alpha(0.5)
            end
            if options.contour_plot
                hold on
                [counts,C] = hist3(data(:,[j i]),[10,10]);
                contour(C{1},C{2},counts')
            end    
            if length(options.true_parameters)>0
                hold on
                % add a marker
                scatter(options.true_parameters(j),options.true_parameters(i),'filled','d')
            end
        end
        set(gca,'fontsize',14) 
       if i==npars
          xlabel(options.parameter_names{j},'Fontsize',15,'fontweight','bold')
       end
       if j==1
           ylabel(options.parameter_names{i},'Fontsize',15,'fontweight','bold')
       end
    end
end
