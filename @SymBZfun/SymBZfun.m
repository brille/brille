% The SymBZfun object holds one 
% py.symbz._symbz.BZGridQ or py.symbz._symbz.BZGridQE object, 
% a function that can be evaluated to fill the object, and one or more
% functions to interpret the interpolation results for, e.g, Horace.
classdef SymBZfun < handle 
    % making copies of the python object is (a) a bad idea for memory 
    % management and (b) possibly bad for memory access --> handle class
    properties
        BZGrid
        filler
        interpreter
    end
    properties (GetAccess = protected, SetAccess = private, SetObservable=true)
        parameterHash = ''
    end
    properties (SetAccess = private)
        isQE = false
        nFill = 1
        shape = {1}
        nRet = 0
        nInt = 1
        rluNeeded = true
        parallel = true
    end
    methods
        function newobj = SymBZfun(BZGrid,varargin)
            kdef = struct('fill',@(x)(1+0*x),...
                          'nfill',[],...
                          'shape',[],...
                          'model','',...
                          'interpret',@(x)(x),...
                          'nret',[],...
                          'rlu',true,...
                          'parallel',true);
            [args,kwds]=symbz.parse_arguments(varargin,kdef,{'rlu'});
            g3type = 'py.symbz._symbz.BZGridQ';  % or py.symbz._symbz.BZGridQcomplex
            g4type = 'py.symbz._symbz.BZGridQE'; % or py.symbz._symbz.BZGridQEcomplex
            if strncmp(class(BZGrid),g4type,length(g4type))
                newobj.isQE=true;
            elseif ~strncmp(class(BZGrid),g3type,length(g3type))
                error('A single %s or %s (or their complex variants) is required as input',g3type,g4type);
            end
            
            fill = kwds.fill;
            if numel(args)>0
                if isa(args{1},'function_handle')
                    fill = args{1};
                elseif iscell(args{1}) && all( cellfun(@(x)(isa(x,'function_handle')), args{1}) )
                    fill = args{1}(1);
                end
            end
            assert(isa(fill,'function_handle'));
            newobj.filler = fill;
            
            % anything that defines 'varargout', including anonymous functions, returns negative nargout
            if ~isempty(kwds.nfill) && isnumeric(kwds.nfill) && isscalar(kwds.nfill)
                nfill = kwds.nfill;
            else
                nfill = abs(nargout(fill)); 
            end
            fshape = kwds.shape; % what is the shape of each filler output
            rlu = kwds.rlu; % does the filler function expect Q in rlu or inverse Angstrom?
            
            nret = [];
            if ~isempty(kwds.model) && ischar(kwds.model)
                switch lower(kwds.model)
                    case 'spinw'
                        nfill = 2;
                        rlu = true;
                        interpret = { @newobj.neutron_spinwave_intensity, @newobj.convolve_modes };
                        nret = [2,1];
                        fshape = {1,[3,3]};
                end
            else
                interpret = kwds.interpret;
                if numel(args)>1
                    interpret = args{2};
                end
            end
            if ~iscell(interpret)
                interpret = {interpret};
            end
            assert( all( cellfun(@(x)(isa(x,'function_handle')), interpret) ),...
                'A single function handle or a cell of function handles is required for the interpreter' );
            newobj.nInt = numel(interpret);
            
            if ~isempty(kwds.nret) && isnumeric(kwds.nret) && numel(kwds.nret)==numel(interpret)
                nret = kwds.nret;
            elseif isempty(nret)
                nret = cellfun(@(x)(abs(nargout(x))),interpret);
            end
            
            if ~iscell(fshape)
                fshape = {fshape};
            end
            assert( ~isempty(fshape) && numel(fshape) == nfill, 'We need to know the shape of the filler output(s)' );
            
            assert( nret(end) == 1, 'the last interpreter function should return a scalar!');
            newobj.nFill = nfill;
            newobj.shape = fshape;
            newobj.rluNeeded = rlu;
            newobj.nRet = nret;
            newobj.interpreter = interpret;
            
            newobj.BZGrid=BZGrid;
            if islogical( kwds.parallel)
                newobj.parallel = kwds.parallel;
            end
        end
        sqw = horace_sqw(obj,qh,qk,ql,en,varargin)
        QorQE = get_mapped(obj)
        fill(obj,varargin)
        intres = interpolate(obj,qh,qk,ql,en)
        sqw = unpolarized_neutron_spinwave_intensity(obj,qh,qk,ql,en,omega,Sab,varargin)
        con = convolve_modes(obj,qh,qk,ql,en,omega,S,varargin)
    end
end
    
