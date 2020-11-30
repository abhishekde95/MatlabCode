function out = sigmoid(in,lambda,offset)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
out = 1./(1 + exp(-lambda*(in-offset)));

end

