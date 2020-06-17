function param_vec = get_param_vec(names,object)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is a function for getting a vector of 
    % parameters from the full parameter vector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    param_vec = [];
    for i=1:length(names)
        param_vec(i) = object.(names{i});
    end
end
