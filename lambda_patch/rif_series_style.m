function [col, ls, mk] = rif_series_style(len, pow)
%RIF_SERIES_STYLE  Distinguishable colors, line styles and markers for the plot series of the patch.
%
%   [col, ls, mk] = rif_series_style(len, pow)
%       len, pow : one entry per series (RIF symbols, power as in the folder name, e.g. 82 for -82 dBm)
%       col      : n x 3 colors. One color per distinct power: up to 10 powers get 10 qualitative colors
%                  (no yellow, readable on white), more powers get an ordered dark-blue -> dark-red map,
%                  so neighbouring powers stay in order in the legend.
%       ls, mk   : n x 1 cells. Line style and marker follow the RIF length:
%                  32 '-' 'o', 64 '--' 's', 128 ':' 'd', 256 '-.' '^' (other lengths '-' 'o').
%   Every (length, power) pair of one analysis therefore looks different, and a condition keeps its look
%   in all figures of that analysis.
%
%   col = rif_series_style(n)   n distinct colors with no length / power meaning (e.g. comb phases).
%
%   MATLAB's lines() has 7 colors and repeats from the 8th series on; the original h0 figure colored by
%   RIF length only, so conditions of one length were indistinguishable.

if nargin == 1 && isscalar(len)
    col = palette(len);
    ls = repmat({'-'}, len, 1); mk = repmat({'o'}, len, 1);
    return;
end
len = len(:); pow = pow(:);
n = numel(len);
pset = unique(pow);                                   % ascending folder power = descending dBm
pc = palette(numel(pset));
col = zeros(n, 3); ls = cell(n, 1); mk = cell(n, 1);
LS = {'-', '--', ':', '-.'}; MK = {'o', 's', 'd', '^'}; LEN = [32 64 128 256];
for i = 1:n
    col(i, :) = pc(pset == pow(i), :);
    q = find(LEN == len(i), 1);
    if isempty(q), q = 1; end
    ls{i} = LS{q}; mk{i} = MK{q};
end
end


% ======================================================================
function c = palette(n)
% blue, orange, green, red, purple, brown, pink, cyan, olive, gray
q = [42 120 214; 235 104 52; 27 150 90; 208 59 59; 130 90 190; ...
     140 86 60; 214 90 170; 20 160 190; 150 150 30; 105 105 105] / 255;
if n <= size(q, 1)
    c = q(1:n, :);
    return;
end
a = [ 30  60 150;  42 120 214;  20 160 160;  60 160  60; 235 140  30; 200  40  40; 120  20  40] / 255;
t = linspace(0, 1, size(a, 1));
c = interp1(t, a, linspace(0, 1, n));
end
