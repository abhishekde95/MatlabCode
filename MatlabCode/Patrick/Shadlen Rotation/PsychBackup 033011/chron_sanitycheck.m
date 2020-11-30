%% This script is to test the working order of 'chronometric.m'
clear all
close all

% Variables
fakedata.coh=[0 .005 .01 .02 .03 .04 .05 .06 .07]'
A=10
B=100
k=-10
th_guess=[1 10 -1]


%% Generate fake data


fakedata.meanrt=(A./(k.*fakedata.coh)).*tanh(A*k.*fakedata.coh)+B;
fakedata.meanrt(isnan(fakedata.meanrt))=A^2+B;
fakedata.meanrt=fakedata.meanrt'
fakedata.sert=.1.*fakedata.meanrt


opts = optimset('MaxFunEvals',10.^6,'MaxIter',10.^6);
[th_fit_chron, prodLt]=fminsearch('chronometric',th_guess,opts,fakedata)
