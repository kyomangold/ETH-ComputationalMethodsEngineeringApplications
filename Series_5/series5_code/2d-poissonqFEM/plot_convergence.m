basename = 'Lshape';
% or e.g.
% basename = 'square'

norm = '';
% or for H1 norm,
% norm = 'H1'

errors = load([basename, '_errors', norm, '.txt']);
numberOfElements = load([basename, '_resolutions.txt']);
diff(log(errors)) ./ diff (log(numberOfElements))
polyfit(log(numberOfElements), log(errors), 1)

loglog(numberOfElements, errors, '-o')

xlabel('Number of free vertices')
ylabel('||u_{FEM}-u||_{L^2}')
grid on
