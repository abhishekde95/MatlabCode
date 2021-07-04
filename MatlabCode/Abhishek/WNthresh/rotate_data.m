function [X_mod,Y_mod] = rotate_data(data_x,data_y,angle)
%Rotate the X and Y data by angle radians
X_mod = data_x*cosd(angle)-data_y*sind(angle);
Y_mod = data_y*cosd(angle)+data_x*sind(angle);


end

