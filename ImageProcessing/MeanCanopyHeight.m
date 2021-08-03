%calculate mean canopy height
clear
load('D:\Mekong_W2015\DataAnalysis\DataReports\Vegetation\2015\Q4_VTAssample.mat')
v = length(vegdat.z);
h = 0.005;
Hmax = h*v; %doesn't include offset from bed of first height slice
mid = h+h/2:h:Hmax+h;
N = vegdat.n;
%calculate h_mean
diff = zeros(v,1);
for i = 1:v-1
    df = N(i)-(N(i+1));
    midpts = mid(i);
    if df < 0
        H(i) = NaN;
        diff(i) = NaN;
    else
        if df == 1
            H(i,1) = midpts;
            diff(i) = df;
        else
            H(i,1:df) = midpts;
            diff(i) = df;
        end
    end
end
H(end,1) = Hmax-h; %include max canopy height
Hz = H(H~=0);Hn = Hz(~isnan(Hz));
H = sort(Hn);
z = h:h:Hmax;
disp(['Mean Canopy height ' num2str(mean(H)) ' m'])
disp(['Stdv. of canopy height ' num2str(std(H)) ' m'])