function bzm = BZMeshQ(BrillouinZone,varargin)
% Create a py.symbz.BZMeshQ object from a py.symbz.BrillouinZone object
% and an optional list of keyword-value paired arguments.
% One argument, 'complex' is used by the MATLAB function to choose a real or
% complex-valued data structure in the mesh; the remaining provided
% keyword-value paired agruments are passed to the Python function.
% As of this writing, the valid keywords are:
%
%   max_size        The maximum edge length of any tetrahedron in the mesh
%                   expressed in inverse Angstrom. Default value -1.0, which
%                   indicates no maximum size.
%   min_angle       The minimum tetrahedron face angle in degrees. Default
%                   value 20.0; -1.0 indicates no minimum. Care should be taken
%                   when setting this value higher than 20 degrees; at some
%                   point the number of mesh points required to satisfy a
%                   large minimum face-angle grows rapidly and above ~30 degrees
%                   there might not be a solution.
%   max_angle       The maximum tetrahedron dihedral angle in degrees. Default
%                   value -1.0, which indicates that an interal default of 179
%                   or 179.99 degrees should be used.
%   min_ratio       The minimum tetrahedron length(?) to circumscribed sphere
%                   radius ratio. Default value -1.0, indicating no limit.
%   max_points      The maximum number of Steiner points to add to the mesh
%                   while attempting to satisfy other provided mesh quality
%                   criteria. This keyword value must be a Python integer and
%                   its default is -1, meaning there is no limit to extra points.
kdef = struct('complex',false);
[args,kwds]=symbz.parse_arguments(varargin,kdef,{'complex'});

reqInType = 'py.symbz._symbz.BrillouinZone';
assert(isa(BrillouinZone,reqInType), ['A single',reqInType,' is required as input']);

if numel(args)>1
    args = pyargs(args{:});
else
    args = pyargs();
end
if kwds.complex
    bzm = py.symbz.BZMeshQcomplex(BrillouinZone, args);
else
    bzm = py.symbz.BZMeshQ(BrillouinZone, args);
end

end
