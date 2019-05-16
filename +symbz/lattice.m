function [dlat,rlat] = lattice(varargin)
% use varargin to define angle units, and whether the parameters describe
% a direct or reciprocal lattice.
kdef = struct('degree',true,'radian',false,'direct',true,'reciprocal',false,'spgr','P 1');
[args,kwds]=symbz.parse_arguments(varargin,kdef,{'degree','radian','direct','reciprocal'});

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

pylens =symbz.m2p( lens(1:3) );
pyangs =symbz.m2p( angs(1:3) );

if kwds.direct && ~kwds.reciprocal
    dlat = py.symbz.Direct(pylens, pyangs, kwds.spgr);
    rlat = dlat.star();
elseif kwds.reciprocal && ~kwds.direct
    rlat = py.symbz.Reciprocal(pylens, pyangs, kwds.spgr);
    dlat = rlat.star();
end
end