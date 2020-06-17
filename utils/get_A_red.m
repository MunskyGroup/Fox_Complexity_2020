function [A_tv,A_constant,k21,A,A4] = get_A_red(parameters,experiment,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a function to get the A matrix for the
% reduced model, where the 3rd state is removed 
% and replaced with bursting dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msize = 3*N;
% state transition matrix
B1 = [-parameters.k12           0                                                  0     ;
       parameters.k12     -parameters.k24-parameters.k22                           0   ; 
       0                        0                                           -parameters.k42-parameters.k44 ]; 

% Get the state transistions.
A1 = kron(eye(N),B1);

% Get the time varying
k21 = @(t) max([0,parameters.k21a-parameters.k21b*experiment.input(t)]);
d1 = repmat([0,1,0],1,N);
A_k21 = spdiags([d1',-d1'],[1,0],msize,msize);

% normal production of RNA
B2 = [parameters.kr1,parameters.kr2, parameters.kr4]; 
d2 = repmat(B2,1,N);
A2 = spdiags([d2',-d2'],[-3,0],msize,msize);

% degradation of RNA 
d3 = repmat(0:N,3)*parameters.g;
A3 = spdiags([d3(:),-d3(:)],[3,0],msize,msize);

% bursting production of RNA
gv = repmat(parameters.p.^(0:N-1)*(1-parameters.p),3,1);
geo_vec = repmat([0,1,0],1,N);
geo_vec_2 = repmat([0,0,1],1,N);
geo_vec = (geo_vec.*gv(:)')';
geo_vec_2 = (geo_vec_2.*gv(:)')';

A4 = zeros(msize,msize);

% There are four ways that bursts are made. k22, k24, k42, k44. 
for i=1:3:N
    A4(i:end,i+1) = parameters.k22*geo_vec(1:end-i+1) + parameters.k24*geo_vec_2(1:end-i+1);
    A4(i:end,i+2) = parameters.k44*geo_vec_2(1:end-i+1)+parameters.k42*geo_vec(1:end-i+1);
end
%Atmp = -1*diag(sum(A4,1));
%A4 = Atmp+A4;

A4(A4<1e-5) = 0; 

A_tv = A_k21;
A_constant = A1+A2+A3+A4;

A = @(t) A1+A2+A3+A_k21*k21(t)+A4;

end
