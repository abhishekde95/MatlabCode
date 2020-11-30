function [val, zmean, zstderr] = hash(x, y);
% function[val, zmean, zstderr] = hash(x,y) outputs val = x#x, zmean = y#x,
% and zstderr = standard error of the set of y values corresponding the the
% same x value

[x1, idx] = sort(x); % x1 = data from array x sorted in ascending order, idx = array of original index position of the now sorted x values
y1 = y(idx); %  y array values resorted to correspond to the ascending x array values now in x1

i = [ find(x1(1:end-1) ~= x1(2:end)) length(x1) ];  % i = array of divider locations between groups of identical values
len = diff([ 0 i ]); % array showing how many elements are numerically identical (each column is a different numerical value) 
val = x1(i); % sorted array of the different numerical values that appear in x, or x1#x1
 
idxlow = [1 cumsum(len(1:end-1))+1]; % array gives the exact indexes where a new numerical value of x element starts
idxhigh = idxlow+len-1; % array gives the exact indexes where the last identical x element value ends
 
m = length(idxlow); % loop will run as many times as there are different x values
for i = 1:m
    zmean(i) = mean(y1(idxlow(i):idxhigh(i))); % gives an array of means of y1 values corresponding to a selected group of identical x1 values, or y#x
    zstderr(i) = std(y1(idxlow(i):idxhigh(i)))/sqrt(len(i));  % gives an array of standard errors of y1 values corresponding to a selected group of idenical x1 values
end