function dockreset()
%DOCKRESET - Close all figures and work-around placement bug.
%
%   When there is a grid of figures and 'close all' is called,
%   a new figure() is placed figure in lower right of grid instead
%   of the upper left.
%
%   This work-around is needed because there are no functions for
%   controlling the placement of figures in a dock window with a grid.
% 
%   Here we activate the first figure using figure(1), pause
%   and then issue 'close all'. The pause is needed for this to 
%   work more reliably (figure execution happens in a different thread
%   and this gives time for the figure(1) command to complete).
%

dock on;
figure(1); 
pause(0.1);
close all;
