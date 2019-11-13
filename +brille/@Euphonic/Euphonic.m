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

classdef Euphonic < handle
  properties
    pyobj
  end
  methods
    function obj = Euphonic(euphonic, varargin)
      kdef = struct('scattering_lengths',[], ...
                    'parallel',true, ...
                    'halfN',[], ...
                    'step',[], ...
                    'unit','rlu', ...
                    'asr',[], ...
                    'precondition',[], ...
                    'set_attrs',[], ...
                    'dipole',[], ...
                    'eta_scale',[], ...
                    'splitting',[], ...
                    'use_primitive',false ...
                    );
      [~,kwds] = brille.parse_arguments(varargin,kdef);
      % Verify that euphonic is a InterpolationData object:
      assert( isa(euphonic,...
          'py.euphonic.data.interpolation.InterpolationData') )
      % Find the key names which have non-empty values:
      keys = fieldnames(kwds);
      i = 1;
      while (i<numel(keys))
          if isempty(kwds.(keys{i}))
              keys(i)=[];
          else
              i = i+1;
          end
      end
      % And store them
      pykwds = cell(2*numel(keys),1);
      for i=1:numel(keys)
        pykwds{2*(i-1)+1} = keys{i};
        pykwds{2*(i-1)+2} = brille.m2p(kwds.(keys{i}));
      end
      obj.pyobj = py.brille.euphonic.FibEu(euphonic, pyargs(pykwds{:}));
    end % intializer
    sqw = horace_sqw(obj,qh,qk,ql,en,varargin)
    wq  = w_q(obj,qh,qk,ql,varargin)
  end % methods
end % classdef
