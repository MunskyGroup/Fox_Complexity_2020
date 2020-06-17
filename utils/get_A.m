function A = get_A(parameters,experiment,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a function to get the A matrix for the
% experiment dsign for the simple 4-state hog model. 
% The A matrix is a function of time. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msize = 4*N;
% state transition matrix
B1 = [-parameters.k12       0                   0                            0;
      parameters.k12   -parameters.k23   parameters.k32                      0; 
       0               parameters.k23   -parameters.k32-parameters.k34   parameters.k43; 
       0                    0               parameters.k34              -parameters.k43];
A1 = kron(eye(N),B1);

% Get the time varying
k21 = @(t) max([0,parameters.k21a-parameters.k21b*experiment.input(t)]);
d1 = repmat([0,1,0,0],1,N);
A_k21 = spdiags([d1',-d1'],[1,0],msize,msize);


% production of RNA
B2 = [parameters.kr1,parameters.kr2,parameters.kr3,parameters.kr4]; 
d2 = repmat(B2,1,N);
A2 = spdiags([d2',-d2'],[-4,0],msize,msize);

% degradation of RNA 
d3 = repmat(0:N,4)*parameters.g;
A3 = spdiags([d3(:),-d3(:)],[4,0],msize,msize);

% Add a reflection boundary condition
%A_k21 = A_k21-spdiags(sum(A_k21)',0,msize,msize);
%A2 = A2-spdiags(sum(A2)',0,msize,msize);
%A3 = A3-spdiags(sum(A3)',0,msize,msize);

A = @(t) A1+A2+A3+A_k21*k21(t);

end
