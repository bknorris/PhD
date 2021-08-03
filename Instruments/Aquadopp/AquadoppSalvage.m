% close
% clear
% dirc = 'D:\Projects\Mekong_F2014\Data\Aquadopp\';
% dd = dir(dirc);
% dname = {dd.name};
% dname = dname(3:9);
% magdec = 0.31;
% for i = 1:length(dname)
%     disp(['Loading data in directory ' dname{i}])
%     cd([dirc dname{i}])
%     ver1 = dir('*_ver1.mat');
%     ver2 = dir('*_ver2.mat');
%     ver1 = {ver1.name};
%     ver2 = {ver2.name};
%     fid = fopen([dirc dname{i} '\metadata.txt']);
%     mdata = textscan(fid,'%s%f','delimiter',',');
%     for ii = 1:length(ver1)
%         if strfind(ver1{ii},'LR')
%             continue
%         else
%             disp(['Loading files ' ver1{ii} ' & ' ver2{ii}])
%             %load the vertwo file, copy ALL fields to a dummy data structure
%             %except for the velocity fields. Then use the start/stop times of
%             %the datenum field to crop the velocity fields of the version 1
%             %file.
%             ornt = mdata{1}(ii);ornt = ornt{:};heading = mdata{2}(ii);
%             DATA = struct();
%             vertwo = load(ver2{ii});fn = fieldnames(vertwo.aqdp);
%             for j = 1:length(fn)
%                 if strcmp(fn{j},'wvel_b1') || strcmp(fn{j},'wvel_b2') || strcmp(fn{j},'wvel_b3') || strcmp(fn{j},'vel_b1') || strcmp(fn{j},'vel_b2') || strcmp(fn{j},'vel_b3')
%                     continue
%                 elseif strcmp(fn{j},'nBursts')
%                     DATA.nSamp = vertwo.aqdp.(fn{j});
%                 else
%                     DATA.(fn{j}) = vertwo.aqdp.(fn{j});
%                 end
%             end
%             start  = DATA.datenum(1);stop = DATA.datenum(end);
%             clear vertwo
%             verone = load(ver1{ii});fn = fieldnames(verone.aqdp);
%             for j = 1:length(fn)
%                 if strcmp(fn{j},'vel_b1') || strcmp(fn{j},'vel_b2') || strcmp(fn{j},'vel_b3')
%                     vel = verone.aqdp.(fn{j});
%                     id = find(verone.aqdp.datenum >= start & verone.aqdp.datenum <= stop);
%                     v = vel(id,:);
%                     DATA.(fn{j}) = v;
%                     clear vel id v
%                 else
%                     continue
%                 end
%             end
%             clear verone
%             %go to the unwrapping step. Unwrap velocities burst-by-burst,
%             %bin-by-bin. First check for bad correlations:
%             nlin = 1E2;maxbr = 1E3;maxg = 1E4;
%             spb = DATA.metadata.spb;
%             int = spb:spb:length(DATA.burst);
%             remd = length(DATA.burst);
%             intv = [1 int remd];
%             for j = 1:length(intv)-1 %burst counter
%                 for k = 1:DATA.nCells; %bin counter
%                     v1 = DATA.vel_b1(intv(j):intv(j+1),k);
%                     v2 = DATA.vel_b2(intv(j):intv(j+1),k);
%                     v3 = DATA.vel_b3(intv(j):intv(j+1),k);
%                     
%                     if strcmp(ornt,'UP') == 1
%                         v1=fliplr(v1);
%                         v2=fliplr(v2);
%                         v3=fliplr(v3);
%                     end
%                     
%                     %Unwrap using Julia's code
%                     vwrap=(max(v1(:))-min(v1(:)))*0.5;
%                     w1=unwrap_w_prof1(v1,vwrap);
%                     vwrap=(max(v2(:))-min(v2(:)))*0.5;
%                     w2=unwrap_w_prof1(v2,vwrap);
%                     vwrap=(max(v3(:))-min(v3(:)))*0.5;
%                     w3=unwrap_w_prof1(v3,vwrap);
%                     
%                     if strcmp(ornt,'UP') == 1
%                         w1=fliplr(w1);
%                         w2=fliplr(w2);
%                         w3=fliplr(w3);
%                     end
%                     
%                     ccrit = 70;
%                     bc1 = find(DATA.cor1(intv(j):intv(j+1),k)<=ccrit);
%                     bc2 = find(DATA.cor2(intv(j):intv(j+1),k)<=ccrit);
%                     bc3 = find(DATA.cor3(intv(j):intv(j+1),k)<=ccrit);
%                     w1(bc1) = NaN;w2(bc2) = NaN;w3(bc3) = NaN;
%                     w1 = cmgbridge(w1,nlin,maxbr,maxg);
%                     w2 = cmgbridge(w2,nlin,maxbr,maxg);
%                     w3 = cmgbridge(w3,nlin,maxbr,maxg);
%                     DATA.wvel_b1(intv(j):intv(j+1),k) = w1;
%                     DATA.wvel_b2(intv(j):intv(j+1),k) = w2;
%                     DATA.wvel_b3(intv(j):intv(j+1),k) = w3;
%                     clear w1 w2 w3 v1 v2 v3 bc1 bc2 bc3
%                 end
%             end
%             
%             disp('Unwrapping complete')
%             tic
%             %Rotate velocities to ENU, then to true N
%             [u,v,w] = aqdpbeam2enu(DATA.nCells,DATA.wvel_b1,DATA.wvel_b2,DATA.wvel_b3,heading,DATA.pitch,DATA.roll,ornt);
%             [u,v] = mag2truenorth(u,v,magdec);
%             DATA.u = u;DATA.v = v;DATA.w = w;
%             disp(['Velocities rotated in ' num2str(toc/60) ' minutes'])
%             
%             %run surface tracking if the instrument is pointing up
%             if strcmp(ornt,'UP')
%                 DATA = aqdpstrack(DATA,DATA.metadata.lat);
%             end
%             
%             %Pad bursts
%             DATA = padbursts(DATA);
%             fields = {'wvel_b1','wvel_b2','wvel_b3'};
%             DATA = rmfield(DATA,fields);
%             
%             aqdp = DATA;
%             clear DATA
%             name = ver1{ii}(1:end-9);
%             save(name,'aqdp','-v7.3')
%             disp(['file ' name '.mat saved'])
%         end
%     end
% end

close
clear
dirc = 'D:\Projects\Mekong_W2015\Data\Aquadopp\';
dd = dir(dirc);
dname = {dd.name};
dname = dname(3:6);
magdec = 0.27;

for i = 4
    disp(['Loading data in directory ' dname{i}])
    cd([dirc dname{i}])
    ver1 = dir('*_ver1.mat');
    ver2 = dir('*_ver2.mat');
    ver1 = {ver1.name};
    ver2 = {ver2.name};
    fid = fopen([dirc dname{i} '\metadata.txt']);
    mdata = textscan(fid,'%s%f','delimiter',',');
    for ii = 5:6
        if strfind(ver1{ii},'LR')
            continue
        else
            disp(['Loading files ' ver1{ii} ' & ' ver2{ii}])
            %load the vertwo file, copy ALL fields to a dummy data structure
            %except for the velocity fields. Then use the start/stop times of
            %the datenum field to crop the velocity fields of the version 1
            %file.
            ornt = mdata{1}(ii);ornt = ornt{:};heading = mdata{2}(ii);
            DATA = struct();
            vertwo = load(ver2{ii});fn = fieldnames(vertwo.aqdp);
            for j = 1:length(fn)
                if strcmp(fn{j},'wvel_b1') || strcmp(fn{j},'wvel_b2') || strcmp(fn{j},'wvel_b3') || strcmp(fn{j},'vel_b1') || strcmp(fn{j},'vel_b2') || strcmp(fn{j},'vel_b3')
                    continue
                elseif strcmp(fn{j},'nBursts')
                    DATA.nSamp = vertwo.aqdp.(fn{j});
                else
                    DATA.(fn{j}) = vertwo.aqdp.(fn{j});
                end
            end
            start  = DATA.datenum(1);stop = DATA.datenum(end);
            clear vertwo
            verone = load(ver1{ii});fn = fieldnames(verone.aqdp);
            for j = 1:length(fn)
                if strcmp(fn{j},'vel_b1') || strcmp(fn{j},'vel_b2') || strcmp(fn{j},'vel_b3')
                    vel = verone.aqdp.(fn{j});
                    id = find(verone.aqdp.datenum >= start & verone.aqdp.datenum <= stop);
                    v = vel(id,:);
                    DATA.(fn{j}) = v;
                    clear vel id v
                else
                    continue
                end
            end
            clear verone
            %go to the unwrapping step. Unwrap velocities burst-by-burst,
            %bin-by-bin. First check for bad correlations:
            nlin = 1E2;maxbr = 1E3;maxg = 1E4;
            spb = DATA.metadata.spb;
            int = spb:spb:length(DATA.burst);
            remd = length(DATA.burst);
            intv = [1 int remd];
            for j = 1:length(intv)-1 %burst counter
                for k = 1:DATA.nCells; %bin counter
                    v1 = DATA.vel_b1(intv(j):intv(j+1),k);
                    v2 = DATA.vel_b2(intv(j):intv(j+1),k);
                    v3 = DATA.vel_b3(intv(j):intv(j+1),k);
                    
                    if strcmp(ornt,'UP') == 1
                        v1=fliplr(v1);
                        v2=fliplr(v2);
                        v3=fliplr(v3);
                    end
                    
                    %Unwrap using Julia's code
                    vwrap=(max(v1(:))-min(v1(:)))*0.5;
                    w1=unwrap_w_prof1(v1,vwrap);
                    vwrap=(max(v2(:))-min(v2(:)))*0.5;
                    w2=unwrap_w_prof1(v2,vwrap);
                    vwrap=(max(v3(:))-min(v3(:)))*0.5;
                    w3=unwrap_w_prof1(v3,vwrap);
                    
                    if strcmp(ornt,'UP') == 1
                        w1=fliplr(w1);
                        w2=fliplr(w2);
                        w3=fliplr(w3);
                    end
                    
                    ccrit = 70;
                    bc1 = find(DATA.cor1(intv(j):intv(j+1),k)<=ccrit);
                    bc2 = find(DATA.cor2(intv(j):intv(j+1),k)<=ccrit);
                    bc3 = find(DATA.cor3(intv(j):intv(j+1),k)<=ccrit);
                    w1(bc1) = NaN;w2(bc2) = NaN;w3(bc3) = NaN;
                    w1 = cmgbridge(w1,nlin,maxbr,maxg);
                    w2 = cmgbridge(w2,nlin,maxbr,maxg);
                    w3 = cmgbridge(w3,nlin,maxbr,maxg);
                    DATA.wvel_b1(intv(j):intv(j+1),k) = w1;
                    DATA.wvel_b2(intv(j):intv(j+1),k) = w2;
                    DATA.wvel_b3(intv(j):intv(j+1),k) = w3;
                    clear w1 w2 w3 v1 v2 v3 bc1 bc2 bc3
                end
            end
            
            disp('Unwrapping complete')
            tic
            %Rotate velocities to ENU, then to true N
            [u,v,w] = aqdpbeam2enu(DATA.nCells,DATA.wvel_b1,DATA.wvel_b2,DATA.wvel_b3,heading,DATA.pitch,DATA.roll,ornt);
            [u,v] = mag2truenorth(u,v,magdec);
            DATA.u = u;DATA.v = v;DATA.w = w;
            disp(['Velocities rotated in ' num2str(toc/60) ' minutes'])
            
            %run surface tracking if the instrument is pointing up
            if strcmp(ornt,'UP')
                DATA = aqdpstrack(DATA,DATA.metadata.lat);
            end
            
            %Pad bursts
            DATA = padbursts(DATA);
            fields = {'wvel_b1','wvel_b2','wvel_b3'};
            DATA = rmfield(DATA,fields);
            
            aqdp = DATA;
            clear DATA
            name = ver1{ii}(1:end-9);
            save(name,'aqdp','-v7.3')
            disp(['file ' name '.mat saved'])
        end
    end
end

% p = aqdp.pressure;
% u = aqdp.u;
% v = aqdp.v;
% w = aqdp.w;
% time = aqdp.datenum;
%
% ax(1) = subplot(211);
% plot(time,p,'k')
% title('Pressure Signal')
% datetickzoom('x','dd HH:MM','keepticks','keeplimits')
% ax(2) = subplot(212);
% plot(time,u,'b'), hold on
% plot(time,v,'r')
% refline(0)
% linkaxes(ax,'x')
% title('Velocities: U - blue, V - red')
% datetickzoom('x','dd HH:MM','keepticks','keeplimits')
%
% ax(1) = subplot(211);
% imagesc(DATA.datenum,DATA.rangebins',DATA.vel_b1')
% caxis([-0.3 0.3])
% datetickzoom('x','dd HH:MM','keepticks','keeplimits')
% title('Without Unwrapping')
% ax(2) = subplot(212);
% imagesc(DATA.datenum,DATA.rangebins',DATA.wvel_b1')
% caxis([-0.3 0.3])
% datetickzoom('x','dd HH:MM','keepticks','keeplimits')
% title('With Unwrapping')
% linkaxes(ax,'xy')