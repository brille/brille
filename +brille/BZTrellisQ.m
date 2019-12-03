% Copyright 2019 Greg Tucker
%
% This file is part of brille.
%
% brille is free software: you can redistribute it and/or modify it under the
% terms of the GNU Affero General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% brille is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.
%
% See the GNU Affero General Public License for more details.
% You should have received a copy of the GNU Affero General Public License
% along with brille. If not, see <https://www.gnu.org/licenses/>.

function bzt = BZTrellisQ(BrillouinZone,varargin)
kdef = struct('max_volume',0.1,'complex',false);
[args,kwds]=brille.parse_arguments(varargin,kdef,{'complex'});

reqInType = 'py.brille._brille.BrillouinZone';
assert(isa(BrillouinZone,reqInType), ['A single',reqInType,' is required as input']);

if kwds.complex
  bzt = py.brille.BZTrellisQcomplex(BrillouinZone, kwds.max_volume);
else
  bzt = py.brille.BZTrellisQ(BrillouinZone, kwds.max_volume);
end

end
