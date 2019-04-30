function bzg = BZGridQE(BrillouinZone,Espec,varargin)
kdef = struct('N',[1,1,1],'d',[0,0,0],'isrlu',true);
[args,kwds]=parse_arguments(varargin,kdef,{'isrlu'});

reqInType = 'py.symbz._symbz.BrillouinZone';
assert(isa(BrillouinZone,reqInType), ['A single',reqInType,' is required as input']);

assert(ismatrix(Espec))
if isa(Espec,'py.numpy.ndarray')
    assert( Espec.size == 3 && Espec.shape(1) == 3 );
else
    assert( numel(Espec) == 3 );
    if Espec(3) < Espec(1)
        tmp = Espec(3);
        Espec(3)=Espec(1);
        Espec(1)=tmp;
    end
    Espec = py.numpy.array( double(Espec(:)) );
end

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
    bzg = py.symbz.BZGridQE(BrillouinZone, Espec, py.numpy.array( double(d(1:3))), isrlu);
else
    bzg = py.symbz.BZGridQE(BrillouinZone, Espec, py.numpy.array( int32(N(1:3))) );
end

end