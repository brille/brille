function bz = brillouinzone(Reciprocal,varargin)
kdef = struct('extent',1);
[args,kwds]=parse_arguments(varargin,kdef);

reqInType = 'py.symbz._symbz.Reciprocal';
assert(isa(Reciprocal,reqInType), ['A single',reqInType,' lattice is required as input']);
if numel(args)>0 && isscalar(args{1})
    extent = int32(args{1});
else
    extent = int32(kwds.extent);
end

bz = py.symbz.BrillouinZone( Reciprocal, extent);
end
