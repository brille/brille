function bzg = BZGridQ(BrillouinZone,varargin)
kdef = struct('N',[1,1,1],'d',[0,0,0],'isrlu',true);
[args,kwds]=symbz.parse_arguments(varargin,kdef,{'isrlu'});

reqInType = 'py.symbz._symbz.BrillouinZone';
assert(isa(BrillouinZone,reqInType), ['A single',reqInType,' is required as input']);

N = kwds.N;
d = kwds.d;
isrlu=kwds.isrlu;
if numel(args)>0
    if ismatrix(args{1}) && numel(args{1})==3
        if sum(mod(args{1},1))==0
            N = args{1};
        else
            d = args{1};
        end
    end
    if isscalar(args{end}) && islogical(args{end})
        isrlu = args{end};
    end
end
if all(d>0)
    bzg = py.symbz.BZGridQ(BrillouinZone, py.numpy.array( double(d(1:3))), isrlu);
else
    bzg = py.symbz.BZGridQ(BrillouinZone, py.numpy.array( uint64(N(1:3))) );
end

end