clear
depname = 'VTA_2_vp3';
savedir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper1\Eps_Vel_Spectra\';
load('D:\Projects\Mekong_W2015\Data\Vectrino\VPTmatrix.mat')
load('d:\Projects\Mekong_W2015\Data\Vector\DPS2\V5109_130315.mat')
% load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\Spectra\VTA_AD5116wvstats.mat')
% plot(wvstats.time,wvstats.hrmsp),datetickzoom('x','dd HH:MM:SS')
% disp(datestr(dat.vpro3.time(end))), %pause
start = datenum(2015,03,14,08,30,00);stop = datenum(2015,03,14,09,00,00);

%calculate principal wave axis with PCA
u = ADV.U;v = ADV.V;
vid = find(ADV.datetime >= start & ADV.datetime <= stop);
u = cmgbridge(u(vid),100,100,1000);v = cmgbridge(v(vid),100,100,1000);
cmp = pca(u,v,[],0);heading = cmp.mdir;
clear ADV

load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\VTA_2vp3Vels.mat')
fn = fieldnames(dat);
for i = 1:length(fn)
    dfn = fieldnames(dat.(fn{i}));
    disp(fn{i})
    id = find(dat.(fn{i}).time >= start & dat.(fn{i}).time <= stop);
    fields = [1:5 8:15];
    for j = 1:length(fields)
        dat.(fn{i}).(dfn{fields(j)}) = dat.(fn{i}).(dfn{fields(j)})(id,:);
    end
    %Transform coordinates in BEAM to XYZ, rotate xyz to principal wave
    %axis from ancillary velocimeter, rotate XYZ back to BEAM.
    T = Tmat.(fn{i});
    b1 = dat.(fn{i}).beam1;b2 = dat.(fn{i}).beam2;
    b3 = dat.(fn{i}).beam3;b4 = dat.(fn{i}).beam4;
    [x,y,z1,z2] = VPTransform(b1,b2,b3,b4,T,'bx');
    %rotate to cross-shore and along-shore
    rot = heading*pi/180;
    x = x.*(ones(size(x))*cos(rot)) + ...
        y.*(ones(size(y))*sin(rot));
    y = -y.*(ones(size(y))*sin(rot)) + ...
        x.*(ones(size(x))*cos(rot));
    [b1,b2,b3,b4] = VPTransform(x,y,z1,z2,T,'xb');
    dat.(fn{i}).beam1 = b1;dat.(fn{i}).beam2 = b2;
    dat.(fn{i}).beam3 = b3;dat.(fn{i}).beam4 = b4;
    %delete XYZ fields to avoid confusion
    dat.(fn{i}) = rmfield(dat.(fn{i}),'x');
    dat.(fn{i}) = rmfield(dat.(fn{i}),'y');
    dat.(fn{i}) = rmfield(dat.(fn{i}),'z1');
    dat.(fn{i}) = rmfield(dat.(fn{i}),'z2');
    dat.(fn{i}).cmt = {'BEAM vels transformed to XYZ,'; 'rotated to PCA,'; 'and retransformed to BEAM.'};
end
dd = datestr(start,'ddHHMMSS');
name = [depname '_' dd '_30min'];
save([savedir name],'-v7.3','dat')