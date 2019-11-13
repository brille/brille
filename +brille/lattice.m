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

function [dlat,rlat] = lattice(varargin)
% use varargin to define angle units, and whether the parameters describe
% a direct or reciprocal lattice.
kdef = struct('degree',true,'radian',false,'direct',true,'reciprocal',false,'spgr','P 1');
[args,kwds]=brille.parse_arguments(varargin,kdef,{'degree','radian','direct','reciprocal'});

assert( numel(args)>0 ,'At least the lattice vector lengths are required to define a lattice.');
lens = args{1};
if numel(args)>1
    angs = args{2};
elseif kwds.degree && ~kwds.radian
    angs = [90,90,90];
elseif kwds.radian && ~kwds.degree
    angs = [1,1,1]*pi/2;
end

assert(numel(lens)>=3 && numel(angs)>=3)
if kwds.degree && ~kwds.radian
    angs = angs / 180 * pi;
end

pylens =brille.m2p( lens(1:3) );
pyangs =brille.m2p( angs(1:3) );

if kwds.direct && ~kwds.reciprocal
    dlat = py.brille.Direct(pylens, pyangs, kwds.spgr);
    rlat = dlat.star();
elseif kwds.reciprocal && ~kwds.direct
    rlat = py.brille.Reciprocal(pylens, pyangs, kwds.spgr);
    dlat = rlat.star();
end
end
