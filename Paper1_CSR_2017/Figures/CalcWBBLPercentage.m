%step 1: calculate difference between no bottom track and bottom track.
fp1 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\TKE\Vertical\';
dir1 = dir([fp1 '*bdadj.mat']);
dir1 = {dir1.name};
fp2 = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\BDTrack\';
dir2 = dir([fp2 '*.mat']);
dir2 = {dir2.name};
avepercent = zeros(4,1);
for j = 1:4
    load([fp1 dir1{j}])
    load([fp2 dir2{j}])
    dn = fieldnames(Stat);
    fn = fieldnames(data);
    vph = [0.062 0.063 0.061];
    percentwbbl = zeros(3,1);
    for i = 1:3
        E1 = Stat.(dn{i}).z1.E;
        vh = vph(i)-0.04-linspace(0.001,0.03,30);
        bd = data.(fn{i}).bd;
        
        %plot
%         figure(i)
%         imagesc(Stat.(dn{i}).time,vh,E1)
%         hold on
%         set(gca,'ydir','normal')
%         plot(data.(fn{i}).time,vph(i)-bd,'r')
        
        %find difference
        numbins = zeros(length(bd),1);
        for ii = 1:length(bd)
            id = find(vph(i)-bd(ii,:) <= vh,1,'last');
            id2 = find(~isnan(E1(:,1)),1,'last');
            if isempty(id) || isempty(id2)
                continue
            else
                if id2 > id
                    numbins(ii) = 0;
                else
                    numbins(ii) = id-id2;
                end
            end
        end
        wbbl = sum(numbins);
        tkedat = numel(~isnan(E1));
        percentwbbl(i) = wbbl/(tkedat+wbbl)*100;
    end
    avepercent(j) = mean(percentwbbl);
end