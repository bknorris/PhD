function rm = runningmean(x,n)
if ~mod(n,2)
    n=n-1;
end

inds=(-(n-1)/2:(n-1)/2)'*ones(1,size(x,1))+ones(n,1)*(1:size(x,1));
inds(find(inds<1))=1;
inds(find(inds>size(x,1)))=size(x,1);

rm=nan*x;
for ii=1:size(x,2)
    thisx=x(:,ii);
    rm(:,ii)=nanmean(thisx(inds))';
end
end
