% plotPSTH.m
%
%	usage: [fig] = PlotPSTH(m, sem, t, varargin)
%	by: Mehrdad Jazayeri
%	last update: 2008-02-16
%    purpose: plots PSTH from 
%
%	varargin: optional args for the plot (e.g. 'k.', or 'Color',[1 1 0])
%
%	varargout:
%
%    example: [fig] = PlotPSTH(m, sem, t, varargin)
%
function [fig] = plotPSTH(m, sem, t, varargin)

%plot(t, m, '.-', varargin{:}); hold on
plot(t, m, '-', varargin{:}); hold on
[errx, erry] = MakeErrorBars(t, m, sem);
plot(errx, erry, '-', varargin{:});
