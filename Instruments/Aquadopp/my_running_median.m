function med=my_running_median(x,n,maxpropnan)

if ~mod(n,2)
    n=n-1;
end

inds=(-(n-1)/2:(n-1)/2)'*ones(1,size(x,1))+ones(n,1)*(1:size(x,1));
inds(find(inds<1))=1;
inds(find(inds>size(x,1)))=size(x,1);

med=nan*x;
if ~exist('maxpropnan')
    for ii=1:size(x,2)
        thisx=x(:,ii);
        med(:,ii)=median(thisx(inds))';
    end
else
    for ii=1:size(x,2)
        thisx=x(:,ii);
        [trsh,propnan]=nanmedian(thisx(inds));
        med(:,ii)=trsh';
        med(find(propnan>maxpropnan),ii)=nan;
    end
end