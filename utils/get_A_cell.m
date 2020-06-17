function [k21,Acell] = get_A_cell(parameters,experiment,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is a function to get the A matrix for the
%experiment dsign for the simple 4-state hog model. 
%The A matrix is a function of time. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msize = 4*N;

% Get A_k12
Bk12 = zeros(4);
Bk12(1,1) = -1;
Bk12(2,1) = 1; 
Acell{1} = sparse(kron(eye(N),Bk12));

% Get A_k23
Bk23 = zeros(4);
Bk23(2,2) = -1;
Bk23(3,2) = 1; 
Acell{3} = sparse(kron(eye(N),Bk23));

% Get A_k32
Bk32 = zeros(4);
Bk32(2,3) = 1;
Bk32(3,3) = -1; 
Acell{4} = sparse(kron(eye(N),Bk32));

% Get A_k34
Bk34 = zeros(4);
Bk34(4,3) = 1;
Bk34(3,3) = -1; 
Acell{5} = sparse(kron(eye(N),Bk34));

% get A_k43
Bk43 = zeros(4);
Bk43(3,4) = 1;
Bk43(4,4) = -1; 
Acell{6} = sparse(kron(eye(N),Bk43));

% production of RNA
B_kr1 = repmat([1 0 0 0 ],1,N);
B_kr2 = repmat([0 1 0 0 ],1,N);
B_kr3 = repmat([0 0 1 0 ],1,N);
B_kr4 = repmat([0 0 0 1 ],1,N);
Acell{7} = spdiags([B_kr1',-B_kr1'],[-4,0],msize,msize);
Acell{8} = spdiags([B_kr2',-B_kr2'],[-4,0],msize,msize);
Acell{9} = spdiags([B_kr3',-B_kr3'],[-4,0],msize,msize);
Acell{10} = spdiags([B_kr4',-B_kr4'],[-4,0],msize,msize);

%degradation of RNA 
d3 = repmat(0:N,4);
Acell{11} = spdiags([d3(:),-d3(:)],[4,0],msize,msize);

%Get the time varying
k21 = @(t) max([0,parameters.k21a-parameters.k21b*experiment.input(t)]);
d1 = repmat([0,1,0,0],1,N);
Acell{2} = spdiags([d1',-d1'],[1,0],msize,msize);

end
