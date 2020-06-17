function Q = get_Q(t,A,Ai,msize)
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Get the RHS for the ODEs
    % for the sensitivity matrix.
    %%%%%%%%%%%%%%%%%%%%%%%%
    Q = Qi+Qitv*k21(t);
    Q = spalloc(2*msize,2*msize,5*msize);
    Q(1:msize,1:msize) = A(t);
    Q(msize+1:end,1:msize) = Ai;
    Q(msize+1:end,msize+1:end) = A(t);
    
end