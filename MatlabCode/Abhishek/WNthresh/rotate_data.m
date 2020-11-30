function [X_mod,Y_mod] = rotate_data(data_x,data_y,angle)
%Rotate the X and Y data by angle radians
X_mod = data_x*cos(angle)-data_y*sin(angle);
Y_mod = data_y*cos(angle)+data_x*sin(angle);


end

