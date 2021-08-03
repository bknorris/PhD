load('D:\Projects\Mekong_W2015\DataAnalysis\Paper3\RS\Two\F2F2_3_RS10min.mat')
bd = load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\F2F2_1_bdtrace.mat');
fn = fieldnames(bd);
bin = [15 15 1];
rho = 1029; %kg/m^3
for i = 1:3
    uw = RS.f2f2.(fn{i}).uw./100;
    t1 = RS.f2f2.(fn{i}).time;
    t2 = bd.(fn{i}).time;
    uw10 = interp1(t1,uw,t2);
    dat.(fn{i}).time = t2;
    dat.(fn{i}).taub = abs(nanmean(uw10(:,1:20))*rho);
end
sdatdir = 'd:\Projects\Mekong_W2015\DataAnalysis\Paper3\RS\Two\';
save([sdatdir 'F2F2_taub'],'-struct','dat')
bd = load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\BottomTrack\F2F3_1_bdtrace.mat');
fn = fieldnames(bd);
bin = [15 15 1];
rho = 1029; %kg/m^3
for i = 1:3
    uw = RS.f2f3.(fn{i}).uw./100;
    t1 = RS.f2f3.(fn{i}).time;
    t2 = bd.(fn{i}).time;
    uw10 = interp1(t1,uw,t2);
    dat.(fn{i}).time = t2;
    dat.(fn{i}).taub = abs(nanmean(uw10(:,1:20))*rho);
end
save([sdatdir 'F2F3_taub'],'-struct','dat')
