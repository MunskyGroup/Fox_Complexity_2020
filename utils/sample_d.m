function pdat = sample_d( Pf, nsample )
    pdat = zeros( size(Pf) );
    w = Pf;
    w(w<0) = 0;
    jsample = randsample( length(Pf), nsample, true, w ); %keyboard
    [ wbin, jbin ] = hist( jsample, unique(jsample) );
    pdat (jbin) = wbin;

end