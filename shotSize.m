function shot_count = shotSize(cell_count)
% Get the total number of counts in cell array structured data
% DKS
% 
% INPUT
% "cellCount" must be a cell of Nshot X Ndim sized array of counts
%
% OUTPUT
% "shot_count" is an array of number of counts in each shot
%

shot_count=cellfun(@(x)size(x,1),cell_count);