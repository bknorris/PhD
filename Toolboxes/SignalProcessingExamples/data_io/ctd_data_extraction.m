close all; clear;

fid = fopen('S27D1.txt','r');

c = 0;
data_start_line = NaN;

while ~feof(fid)

    c = c + 1;    
    
    curr_line = fgets(fid);

    if (curr_line(1:5)=='*END*')
        data_start_line = c + 1;
    end

    if (c >= data_start_line-1)
        data_tmp = fscanf(fid,'%f');
    end

end

leng = length(data_tmp);
data_records = leng/7;
data = reshape(data_tmp,7,data_records);
data = data';

fclose(fid);
