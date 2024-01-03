basename = 'Lshape';
% or e.g.
% basename = 'square'
vertices = load([basename, '_vertices.txt']);
indices = load([basename, '_triangles.txt']);
indices = indices + 1; % adjust for the 1-indexing in MATLAB
values = load([basename, '_values.txt']);

disp(['Using ', num2str(length(vertices)), ' vertices']);

trimesh(indices, vertices(:,1), vertices(:,2), values)