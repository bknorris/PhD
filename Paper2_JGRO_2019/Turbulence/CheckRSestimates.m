%Check Reynolds stress against small adjustments to the orientation of the
%instrument. RS estimates are purportedly affected by small changes in
%instrment tilt. Determine the accuracy of the RS estimates by comparing
%these to other estimates computed from 5deg pitch/roll adjustments to the
%velocities. 
clear
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\V3\RStress_10min_v3_1deg_pitch.mat')
pitchp = RS;
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\V3\RStress_10min_v3_-1deg_pitch.mat')
pitchn = RS;
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\V3\RStress_10min_v3_1deg_roll.mat')
rollp = RS;
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\V3\RStress_10min_v3_-1deg_roll.mat')
rolln = RS;
% load('D:\Projects\Mekong_W2015\DataAnalysis\Paper1\RStress_10min_-5deg.mat')
% minus5 = RS;
load('D:\Projects\Mekong_W2015\DataAnalysis\Paper2\RS_estimates\V3\RStress_10min_v3.mat')

disp('Listing percent changes from 5% pitch & roll adjustment to velocities');
dn = fieldnames(RS);
fn = fieldnames(RS.day1);
for i = 1:4
    disp(dn{i})
    for j = 1:3
        disp(fn{j})
        uwci = zeros(35,5);
        vwci = zeros(35,5);
        ppuw = pitchp.(dn{i}).(fn{j}).uw./100;
        ppvw = pitchp.(dn{i}).(fn{j}).vw./100;
        pnuw = pitchn.(dn{i}).(fn{j}).uw./100;
        pnvw = pitchn.(dn{i}).(fn{j}).vw./100;
        rpuw = rollp.(dn{i}).(fn{j}).uw./100;
        rpvw = rollp.(dn{i}).(fn{j}).vw./100;
        rnuw = rolln.(dn{i}).(fn{j}).uw./100;
        rnvw = rolln.(dn{i}).(fn{j}).vw./100;
        uw = RS.(dn{i}).(fn{j}).uw./100;
        vw = RS.(dn{i}).(fn{j}).vw./100;
        %UW
        uwci(:,1) = 1.95*(nanstd(ppuw)./sqrt(length(ppuw)));
        uwci(:,2) = 1.95*(nanstd(pnuw)./sqrt(length(pnuw)));
        uwci(:,3) = 1.95*(nanstd(rpuw)./sqrt(length(rpuw)));
        uwci(:,4) = 1.95*(nanstd(rnuw)./sqrt(length(rnuw)));
        uwci(:,5) = 1.95*(nanstd(uw)./sqrt(length(uw)));
        m = max(uwci(:,1:4),[],2);
        totalci = abs(uwci(:,5)-m);
        fprintf('Maximum CI in bins 1-35 for uw stresses: %0.2d Pa\n',nanmax(totalci))
%         fprintf('Mean of bins 1-35 for uw stresses: %0.2d Pa\n',nanmean(nanmean(abs(uw))))
        fprintf('Percentage of total of error estimates to mean values: %0.2f%%\n\n',(nanmax(totalci))./nanmean(nanmean(abs(uw)))*100)

        %VW
%         vwci(:,1) = 1.95*(nanstd(ppvw)./sqrt(length(ppvw)));
%         vwci(:,2) = 1.95*(nanstd(pnvw)./sqrt(length(pnvw)));
%         vwci(:,3) = 1.95*(nanstd(rpvw)./sqrt(length(rpvw)));
%         vwci(:,4) = 1.95*(nanstd(rnvw)./sqrt(length(rnvw)));
%         vwci(:,5) = 1.95*(nanstd(vw)./sqrt(length(vw)));
%         m = max(vwci(:,1:4),[],2);
%         totalci = abs(vwci(:,5)-m);
%         fprintf('Maximum CI in bins 1-35 for vw stresses: %0.2d Pa\n',nanmax(totalci))
% %         fprintf('Mean of bins 1-35 for vw stresses: %0.2d Pa\n',nanmean(nanmean(abs(vw))))
%         fprintf('Percentage of total of error estimates to mean values: %0.2f%%\n\n',(nanmax(totalci))./nanmean(nanmean(abs(vw)))*100)
    end
end
        