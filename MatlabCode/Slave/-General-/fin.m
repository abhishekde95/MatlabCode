save_breakpoints(dbstatus());
close all hidden
clear all
clear classes
dbstop(save_breakpoints());
munlock('save_breakpoints');
clc
