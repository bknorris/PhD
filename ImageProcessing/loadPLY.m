function ply = loadPLY(path,name)
fileName = strcat(path, name);
disp(['Loading ' fileName '.ply']);

idx = 1;
fid = fopen(fileName, 'r');
while isempty(strfind(fgets(fid), 'end_header'))
    idx = idx + 1;
end
fclose(fid);

data = textread(fileName, '%s','delimiter', '\n');
data = data(idx+1:length(data),1);
data = (cellfun(@(x) strread(x,'%s','delimiter',' '), data, 'UniformOutput', false));

if isempty(data{length(data)})
    data(length(data))=[];
end

ply = str2double([data{:}].');
ply = ply(:,1:6);

save(strcat(path,name(1:numel(name)-4),'.mat'), 'ply');

