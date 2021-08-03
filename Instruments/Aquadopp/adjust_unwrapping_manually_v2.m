%%%% if red when should be blue 

n=8; % n/2 = number of 'stripes' to change


q=ginput(n)

for kk=1:length(q)/2
in=round(q(2*kk-1,1));
out=round(q(kk*2,1));    


for ii=in:out
    for jj=1:28
        if aquadopp_burst.vel_b1_uw(ii,jj)>0
            aquadopp_burst.vel_b1_uw(ii,jj)=aquadopp_burst.vel_b1_uw(ii,jj)-2*vwrap;
        end
    end
end

end


%%%% if blue when should be red

n=8; % n/2 = number of 'stripes' to change


q=ginput(n)

for kk=1:length(q)/2
in=round(q(2*kk-1,1));
out=round(q(kk*2,1));    


for ii=in:out
    for jj=1:28
        if aquadopp_burst.vel_b1_uw(ii,jj)<0
            aquadopp_burst.vel_b1_uw(ii,jj)=aquadopp_burst.vel_b1_uw(ii,jj)+2*vwrap;
        end
    end
end

end


% Click whole regions to NaN out 


 
q=ginput(2)

for kk=1:length(q)/2
    in=round(q(2*kk-1,1));
    out=round(q(kk*2,1)); 
    aquadopp_burst.vel_b1_uw(in:out,:)=NaN;
end
 
  
