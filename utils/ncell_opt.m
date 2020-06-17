function [metric_val,FIM] = ncell_opt(x,FIMs,metric)
    FIM = zeros(size(FIMs{1}));
    for k = 1:length(FIMs)
        FIM = FIM+FIMs{k}*x(k);
    end
    switch metric
        case 'degradation' 
            metric_val = FIM(end,end);
        case 'subset_volume'
%             keepers = [1,2,5,6,9];
%             sub_FIM = FIM(keepers,keepers);
            metric_val = prod(eig(FIM));
    end
end
