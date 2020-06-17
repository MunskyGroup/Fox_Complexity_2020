function [marginal,pout] = solve_fsp_tv(model)

    % Get the A matrix
    A = get_A(model.parameters,model.experiment,model.N);
    J = @(t,y) A(t);
    Jsp = A(100)>0;
    % Initial condition
    p0 = zeros(4*model.N,1);
    p0(1) = 1;

    tic;
    options = odeset('Jacobian',J,'JPattern',Jsp); 
    ODE = @(t,x) A(t)*x;
    [tout,pout] = ode23s(ODE,model.tvec,p0,options);
    toc

    % get the marginal distribution of RNA 
    marginal=zeros(length(tout),model.N);
    for i=1:length(tout)     
        marginal(i,:) = sum(reshape(pout(i,:),4,model.N));
    end


end
