clear,close all
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_flood.mat');
data.flood = dat;
load('d:\Projects\Mekong_W2015\DataAnalysis\Paper3\Wavelet\CombEventdata_ebb.mat');
data.ebb = dat;clear dat


bfmd = abs([data.flood.mud.wave.deltbd; data.flood.mud.ig.deltbd]);
id1 = find(bfmd<=0.001);
bffr = abs([data.flood.fringe.wave.deltbd; data.flood.fringe.ig.deltbd]);
id2 = find(bffr<=0.001);
bffo = abs([data.flood.forest.wave.deltbd; data.flood.forest.ig.deltbd]);
id3 = find(bffo<=0.001);

bemd = abs([data.ebb.mud.wave.deltbd; data.ebb.mud.ig.deltbd]);
id4 = find(bemd<=0.001);
befr = abs([data.ebb.fringe.wave.deltbd; data.ebb.fringe.ig.deltbd]);
id5 = find(befr<=0.001);
befo = abs([data.ebb.forest.wave.deltbd; data.ebb.forest.ig.deltbd]);
id6 = find(befo<=0.001);

fprintf('Flood- Percent values within +/-0.001 m, Mudflat: %0.1f%%\n',(numel(id1)/numel(bfmd))*100)
fprintf('Flood- Percent values within +/-0.001 m, Fringe: %0.1f%%\n',(numel(id2)/numel(bffr))*100)
fprintf('Flood- Percent values within +/-0.001 m, Forest: %0.1f%%\n',(numel(id3)/numel(bffo))*100)
fprintf('Ebb- Percent values within +/-0.001 m, Mudflat: %0.1f%%\n',(numel(id4)/numel(bemd))*100)
fprintf('Ebb- Percent values within +/-0.001 m, Fringe: %0.1f%%\n',(numel(id5)/numel(befr))*100)
fprintf('Ebb- Percent values within +/-0.001 m, Forest: %0.1f%%\n',(numel(id6)/numel(befo))*100)
