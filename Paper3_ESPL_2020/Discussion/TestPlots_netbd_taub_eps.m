clear, close all
load('e:\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventData_ebb_v2.mat');
flds = fieldnames(dat);
tcr = [0.18,0.26,0.32];
symb = {'+','*','x'};
cl = flipud([0.1 0.1 0.1;0.5 0.5 0.5;0.7 0.7 0.7]);


f1 = figure(1);
set(f1,'PaperOrientation','portrait',...
    'position',[400 100   500   400],...
    'renderer','painters');
for i = 1:length(flds)
    deltbd = [dat.(flds{i}).wave.deltbd; dat.(flds{i}).ig.deltbd];
    taub = [dat.(flds{i}).wave.taub; dat.(flds{i}).ig.taub];
    eps = [dat.(flds{i}).wave.eps; dat.(flds{i}).ig.eps];
    
    %try finding bed states that are greater than the bed shear stress?
    id = find(taub>tcr(i));
    semilogx(eps(id),taub(id),...
        'marker',symb{i},'color',cl(i,:),...
        'markersize',6,'linestyle','none'),hold on
    
    
%     subplot(121)
%     x = eps(id);y = abs(deltbd(id));
%     id2 = find(abs(deltbd(id))<0.03);
%     semilogx(x(id2),y(id2),...
%         'marker',symb{i},'color',cl(i,:),...
%         'markersize',6,'linestyle','none'),hold on
%     subplot(122)
%     x = taub(id);y = abs(deltbd(id));
%     id2 = find(abs(deltbd(id))<0.03);
%     plot(x(id2),y(id2),...
%         'marker',symb{i},'color',cl(i,:),...
%         'markersize',6,'linestyle','none'),hold on
end


