function maxim=runningmax(x,n,maxpropnan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs: x - dataset
%         n - interval
%         maxpropnan - ??????
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~mod(n,2)
    n=n-1;
end

inds=(-(n-1)/2:(n-1)/2)'*ones(1,size(x,1))+ones(n,1)*(1:size(x,1));
inds(find(inds<1))=1;
inds(find(inds>size(x,1)))=size(x,1);

maxim=nan*x;
if ~exist('maxpropnan')
    for ii=1:size(x,2)
        thisx=x(:,ii);
        maxim(:,ii)=max(thisx(inds))';
    end
else
    for ii=1:size(x,2)
        thisx=x(:,ii);
        [trsh,propnan]=nanmax(thisx(inds));
        maxim(:,ii)=trsh';
        maxim(find(propnan>maxpropnan),ii)=nan;
    end
end