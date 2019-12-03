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

% The Brille object holds one
% py.brille._brille.BZGridQ or py.brille._brille.BZGridQE object,
% a function that can be evaluated to fill the object, and one or more
% functions to interpret the interpolation results for, e.g, Horace.
classdef Brille < handle
    % making copies of the python object is (a) a bad idea for memory
    % management and (b) possibly bad for memory access --> handle class
    properties
        pygrid
        filler
        interpreter
    end
    properties (GetAccess = protected, SetAccess = private, SetObservable=true)
        parameterHash = ''
    end
    properties (SetAccess = private)
        isQE = false
        nFill = 1
        nFillers = 1;
        span  = 1
        shape = {1}
        nRet = 0
        nInt = 1
        rluNeeded = true
        parallel = true
        formfact = false
        magneticion
        formfactfun
        Qscale = eye(4);
        Qtrans = eye(4);
    end
    methods
        function obj = Brille(ingrid,varargin)
            kdef = struct('fill',@(x)(1+0*x),...
                          'nfill',[],...
                          'shape',[],...
                          'model','',...
                          'interpret',@(x)(x),...
                          'nret',[],...
                          'rlu',true,...
                          'parallel',true,...
                          'formfact',false,...
                          'magneticion','',...
                          'formfactfun',@sw_mff,...
                          'Qscale',eye(4),...
                          'Qtrans',eye(4));
            [args,kwds]=brille.parse_arguments(varargin,kdef,{'rlu'});
            g3type = 'py.brille._brille.BZGridQ';  % or py.brille._brille.BZGridQcomplex
            g4type = 'py.brille._brille.BZGridQE'; % or py.brille._brille.BZGridQEcomplex
            m3type = 'py.brille._brille.BZMeshQ';
            n3type = 'py.brille._brille.BZNestQ';
            t3type = 'py.brille._brille.BZTrellisQ';
            if strncmp(class(ingrid),g4type,length(g4type))
                obj.isQE=true;
            elseif strncmp(class(ingrid),g3type,length(g3type))
                warning('The BZGridQ and BZGridQcomplex objects are inefficient. Consider using a BZTrellisQ* instead.');
            elseif strncmp(class(ingrid),m3type,length(m3type))
                warning('The BZMeshQ and BZMeshQcomplex objects are inefficient. Consider using a BZNestQ* instead.');
            elseif strncmp(class(ingrid),n3type,length(n3type))
                warning('The BZNestQ and BZNetstQcomplex objects are slow to locate points. Consider using a BZTrellisQ* instead.');
            elseif ~strncmp(class(ingrid),t3type,length(t3type))
                error('Expected input of a py.brille._brille.BZ{Grid,Mesh,Nest,Trellis}Q* object')
            end
            if islogical( kwds.parallel)
                obj.parallel = kwds.parallel;
            end
            if islogical( kwds.formfact )
                obj.formfact = kwds.formfact;
            end
            if obj.formfact && ~isempty(kwds.magneticion)
                obj.magneticion = kwds.magneticion;
            end
            if obj.formfact && isa(kwds.formfactfun,'function_handle')
                obj.formfactfun = kwds.formfactfun;
            end
            if isnumeric(kwds.Qscale) && ismatrix(kwds.Qscale)
                if numel(kwds.Qscale)==16
                    obj.Qscale(:) = kwds.Qscale(:);
                elseif numel(kwds.Qscale)==9
                    obj.Qscale([1,2,3,5,6,7,9,10,11])=kwds.Qscale(:);
                end
            end
            if isnumeric(kwds.Qtrans) && ismatrix(kwds.Qtrans)
                if numel(kwds.Qtrans)==16
                    obj.Qtrans(:) = kwds.Qtrans(:);
                elseif numel(kwds.Qtrans)==9
                    obj.Qtrans([1,2,3,5,6,7,9,10,11])=kwds.Qtrans(:);
                end
            end

            if iscell(kwds.fill)
                fill = kwds.fill;
            else
                fill = {kwds.fill};
            end
            if numel(args)>0
                if isa(args{1},'function_handle')
                    fill = args(1);
                elseif iscell(args{1}) && all( cellfun(@(x)(isa(x,'function_handle')), args{1}) )
                    fill = args{1};
                end
            end
            assert(iscell(fill) && all( cellfun(@(x)(isa(x,'function_handle')), fill) ));
            obj.filler = fill;
            obj.nFillers = length(fill);

            % anything that defines 'varargout', including anonymous functions, returns negative nargout
            if ~isempty(kwds.nfill) && isnumeric(kwds.nfill) && isscalar(kwds.nfill)
                nfill = kwds.nfill;
            else
                nfill = abs(nargout(fill{1}));
            end
            fshape = kwds.shape; % what is the shape of each filler output
            rlu = kwds.rlu; % does the filler function expect Q in rlu or inverse Angstrom?

            nret = [];
            if ~isempty(kwds.model) && ischar(kwds.model)
                switch lower(kwds.model)
                    case 'spinw'
                        nfill = 2;
                        % if obj.nFillers==1
                        %     obj.filler = [ obj.filler {@brille.modesort} ];
                        %     obj.nFillers=2;
                        % end
                        rlu = true;
                        interpret = { @obj.neutron_spinwave_intensity, @obj.convolve_modes };
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
            obj.nInt = numel(interpret);

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
            obj.nFill = nfill;
            obj.shape = fshape;
            obj.rluNeeded = rlu;
            obj.nRet = nret;
            obj.interpreter = interpret;

            obj.pygrid=ingrid;
        end
        sqw = horace_sqw(obj,qh,qk,ql,en,varargin)
        QorQE = get_mapped(obj)
        fill(obj,varargin)
        intres = interpolate(obj,qh,qk,ql,en)
        sqw = unpolarized_neutron_spinwave_intensity(obj,qh,qk,ql,en,omega,Sab,varargin)
        con = convolve_modes(obj,qh,qk,ql,en,omega,S,varargin)
    end
end
