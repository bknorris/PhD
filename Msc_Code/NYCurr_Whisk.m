function NYCurr_Whisk(u,v,time,data_scale,plot_scale,range_y);
%   NYCurr_Whisk(u,v,time,data_scale,plot_scale,y_range);
%
%Inputs:  u (eastward component of velocity), v (northward component of velocity)
%data_scale (some scaling factor of the data - increasing this number increases
%the size of the whiskers)
%
%Program designed to plot time series current velocity data
%as "whiskers" coming off of the x-axis.  So, for each time step, a line
%is drawn off of the x axis of the proper magnitude and angle.  The y-axis
%is labeled in the units of the velocity and the x-axis is labeled in time
%units
%
%KEEP IN MIND THAT IN ORDER FOR THESE WHISKERS TO BE AT THE PROPER ANGLE
%AND MAGNITUDE, THE AXES MUST BE SET AS EQUAL (THE PLOT DOES THIS).  THIS
%IS WHY 'DATA_SCALE' AND 'PLOT_SCALE' INPUTS ARE NECESSARY:  OTHERWISE, YOU
%WILL GET A PLOT WHICH MAY BE SQUISHED IN THE X OR Y DIRECTION.  THESE ARE
%BUILT INTO THE PROGRAM

range_x = [yearday(2013,12,01) yearday(2015,4,20)];


%How many whiskers are needed
num_whiskers = length(u);

%Scale the data & the time scale
u_scale = data_scale.*u;
v_scale = data_scale.*v;
plot_loc = plot_scale.*time;

hold on

for index_hit = 1:num_whiskers;
plot([plot_loc(index_hit); u_scale(index_hit)+plot_loc(index_hit)],[0; v_scale(index_hit)]);
end

min_x = range_x(1);
max_x = range_x(2);
min_y = range_y(1);
max_y = range_y(2);

%Can be tweaked here to change tick locations
rangeofy = max_y-min_y;
space_of_ticks = floor(rangeofy/4);


axis equal
axis([min_x*plot_scale max_x*plot_scale min_y*data_scale max_y*data_scale])
set(gca,'ytick',[min_y:space_of_ticks:max_y].*data_scale);
set(gca,'yticklabel',[min_y:space_of_ticks:max_y])
l = line([min_x*plot_scale max_x*plot_scale],[0 0]);
set(l,'color','k')
gregaxd_fooled(time,plot_scale,5);
box('on')