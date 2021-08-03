function [v1,v2,v3]=interpolate_bad_correlations(v1,v2,v3,c1,c2,c3,rb,ccrit) 

if ~exist('ccrit')
    ccrit=70;
end


   for ii=1:size(v1,1)
        clear goodind* badind* 
        goodind1=find(c1(ii,:)>=ccrit); % Find high correlation points for each beam
        goodind2=find(c2(ii,:)>=ccrit);
        goodind3=find(c3(ii,:)>=ccrit);

        % Interpolate
        if length(goodind1)>=2 && length(goodind2)>=2 && length(goodind3)>=2 
            v1(ii,:)=interp1(rb(goodind1)...
                ,v1(ii,goodind1),rb,'nearest','extrap');
            v2(ii,:)=interp1(rb(goodind2)...
                ,v2(ii,goodind2),rb,'nearest','extrap');
            v3(ii,:)=interp1(rb(goodind3)...
                ,v3(ii,goodind3),rb,'nearest','extrap');
        end
   end
   
  
   