% Copyright 2019 Greg Tucker
%
% This file is part of fibril.
%
% fibril is free software: you can redistribute it and/or modify it under the
% terms of the GNU Affero General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% fibril is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.
%
% See the GNU Affero General Public License for more details.
% You should have received a copy of the GNU Affero General Public License
% along with fibril. If not, see <https://www.gnu.org/licenses/>.

function bzg = BZGridQ(BrillouinZone,varargin)
kdef = struct('N',[1,1,1],'d',[0,0,0],'isrlu',true,'complex',false);
[args,kwds]=fibril.parse_arguments(varargin,kdef,{'isrlu','complex'});

reqInType = 'py.fibril._fibril.BrillouinZone';
assert(isa(BrillouinZone,reqInType), ['A single',reqInType,' is required as input']);

N = kwds.N;
d = kwds.d;
isrlu=kwds.isrlu;
if numel(args)>0
    if ismatrix(args{1}) && numel(args{1})==3 || (isa(args{1},'py.numpy.ndarray')&&args{1}.size==3)
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
if ~isa(N,'py.numpy.ndarray')
    N = fibril.m2p(N, 'uint64');
end
if ~isa(d,'py.numpy.ndarray')
    d = fibril.m2p(d, 'double');
end
if all(d>0)
    if kwds.complex
        bzg = py.fibril.BZGridQcomplex(BrillouinZone, d, isrlu);
    else
        bzg = py.fibril.BZGridQ(BrillouinZone, d, isrlu);
    end
else
    if kwds.complex
        bzg = py.fibril.BZGridQcomplex(BrillouinZone, N );
    else
        bzg = py.fibril.BZGridQ(BrillouinZone, N );
    end
end

end
