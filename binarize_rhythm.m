function binary = binarize_rhythm(pattern)

% function binaries = binarize_rhythm(pattern)
% 
% Input:
%   pattern = pattern to binarize (vector)
% 
% Output:
%   binary = 0s and 1s for input pattern. Can be used to generate sound or
%   as input to Povel & Essens model, for example.

binary = [];
for ii = 1:length(pattern)
    binary = [binary 1 zeros(1,pattern(1,ii)-1)];
end