classdef SymSim < handle
  properties
    pyobj
  end
  methods
    function obj = SymSim(pysimphony, varargin)
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
                    'splitting',[] ...
                    );
      [~,kwds] = symbz.parse_arguments(varargin,kdef);
      % Verify that pysimphony is a InterpolationData object:
      assert( isa(pysimphony,...
          'py.simphony.data.interpolation.InterpolationData') )
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
        pykwds{2*(i-1)+2} = symbz.m2p(kwds.(keys{i}));
      end
      obj.pyobj = py.symbz.simphony.SymSim(pysimphony, pyargs(pykwds{:}));
    end % intializer
    sqw = horace_sqw(obj,qh,qk,ql,en,varargin)
  end % methods
end % classdef
