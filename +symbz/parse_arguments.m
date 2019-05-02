function [par,keyval,present,filled,ok,mess]=parse_arguments(args,varargin)
% Parse a cell array to find values of parameters and named keywords
%
% Basic use:
% ----------
% Unlimited number of parameters:
%   >> [par,keyval,present,filled,ok,mess]=...
%              symbz.parse_arguments (args, keyval_def)
%
% Specify the number of required parameters and number of optional parameters:
%   >> [...] = symbz.parse_arguments (args, npar_req, npar_opt, keyval_def)
%
% Specify the names of required parameters as a cell array of strings, and
% optional parameters as a structure giving names and values:
%   >> [...] = symbz.parse_arguments (args, par_req, par_opt, keyval_def)
%
%
% Optional additional arguments (one or both, in either order) 
% -----------------------------
% Indicate which keywords are logical flags:
%   >> [...] = symbz.parse_arguments (..., flagnames)
%
% Additional options (contained in a structure):
%   >> [...] = symbz.parse_arguments (..., opt)
%
%
% Input:
% ------
%   args        Cell array of arguments to be parsed. The list of arguments is
%              in general a set of parameters followed by optional keyword-value
%              pairs:
%                 par_1, par_2, par_3, ..., name_1, val_1, name_2, val_2, ...
%
%               where
%                 par_1, par_2, ... are the values of positional parameters
%
%                 name_1 ,val_1, ... are optional named arguments and their
%                                     corresponding values.
%               Note:
%               - The valid keywords are the field names of the structure
%                keyval_def (details below), and their default output values
%                are given by the values of those fields.
%
%               - The keyword names can be exact matches or unambiguous
%                abbreviations [exact only if opt.keys_exact ia set to true].
%
%               - Keyword names that are also defined as logical flags in
%                flagnames (below) can also be given as their negation
%                by default [change this by setting opt.flags_noneg==true]
%                 e.g. if 'mask' is a keyword, 'nomask' is also valid.
%
%               - The value of a logical flag does not need to be given a
%                value: if 'foo' is a flag then
%                   ...,'foo',1,...          sets argout.foo = 1
%                   ...,'foo',...            sets argout.foo = 1
%                   ...,'foo',0,...          sets argout.foo = 0
%                   ...,'nofoo',...          sets argout.foo = 0
%
%               - If a keyword is not given, then the default value defined by
%                keyval_def is returned.
%
%               - A prefix to the keywords can be specified in the options
%                (see below). For example if a field of keyval_def is
%                'mask' and the prefix '-' is specified then the keyowrd
%                must appear in args as '-mask', and its negatin (if it is
%                a logical flag) is '-nomask'.
%
%               - By default, all positional parameters must appear first; the
%                first time a character string is encountered that matches a
%                parameter name that marks the point when only keyword-value
%                pairs and logical flags can follow. The default can be changed
%                to allow mixed parameters and keywords if opt_keys_at_end is
%                set to false.
%
% Optionally: give either the number of required and optional positional
% parameters:
%   npar_req    The number of required positional parameters. It can be 
%              0,1,2,...  [Default: =[] which is equivalent to 0]
%
%   npar_opt    The number of optional positional parameters. It
%              can be 0,1,2,... Inf . [Default: =[] which is equivalent to 0]
%
% Or the alternative option: give the names of required positional parameters
% and both the names and values of optional positional parameters:
%   par_req     Cell array with the names of required positional parameters
%              [Default: =[] which means no required parameters]
%
%   par_opt     Structure with the names and values of optional positional
%              parameters. The names of the parameters are given by the
%              field names of the structure, and the values by the values of
%              the fields. [Default=[] which means no optional parameters]
%
% If neither option is given, then this is interpreted as no required 
% parameters and unlimited optional parameters i.e. it is the same as
% npar_req=0, npar+opt=Inf [Note this is NOT the same as the default
% numeric input]
%
%               
%   keyval_def  Structure with field names giving the parameter names, and
%              default values. Note that some keywords may not be permitted
%              if negation of logical flags is permitted (which is the
%              default behaviour)  e.g. if 'mask' is a logical flag then
%              'nomask' is implicity also a name, and this is not permiitted
%              to be a field of keyval_def.
%
%   flagnames   [Optional] Cell array giving the keywords that can only be
%              logical 0 or 1. The names must be given in full (i.e. not
%              abbreviations), and the corresponding default values given in
%              keyval_def must be false, true, 0 or 1, otherwise an error
%              error is returned.
%               If empty, or omitted, then no arguments will be logical flags.
%
%   opt         [Optional] Structure giving control options. Valid fields are
%              given below, together with the default values if the field is
%              not given.
%               prefix      Value of prefix to keywords. Must be a character
%                          string e.g. '-'. [Default: '' i.e. no prefix]
%
%               prefix_req  True if require that the prefix be present, false
%                          otherwise. [Default: true]
%                           For example if 'mask' is a keyword, opt.prefix='-'
%                          opt.prefix_req=false, then 'mask' and '-mask' are
%                          both valid.
%
%               flags_noneg True if flags cannot be negated e.g. if 'full' is
%                          a flag, 'nofull' is not permitted. [Default: false]
%
%               flags_noval True if flags cannot be given values e.g.
%                          ...,'foo',1,... is not permitted; the value is
%                          set by specifing ...,'foo',... or ...,'nofoo',...
%                           [Default: false]
%                           Note: opt.flags_noneg and opt.flags_noval cannot
%                          both be true if a logical flag has default value
%                          true in keyval_def, because then there is no way
%                          to set its value to false.
%
%               keys_exact  True if exact match to keywords is required
%                          [Default: false]
%
%               keys_at_end True if keywords must appear at the end of the
%                          argument list; otherwise keywords and un-named
%                          parameters can be mixed. [Default: true]
%
%               keys_once   True if keywords are only allowed to appear once
%                          in the argument list [Default: true]
%                           If false i.e. keywords can be repeated, then the
%                          last occurence takes precedence.
%
%               noffset     Offset (>=0) for error message display [Default: 0]
%                          If not all arguments are passed to symbz.parse_arguments
%                          then if an error is found at the third position
%                          in args, the error message that symbz.parse_arguments
%                          gives will be stated at the third position, but
%                          this will not be at the true position in the list.
%                          Give the offset here. For example, you might not
%                          pass the first three arguments because these must
%                          always be present, set opt.noffset=3.
%
% Output:
% -------
% par       If only the number of required and optional positional parameters
%          was given (or neither option): Cell array that contains values of
%          arguments that do not correspond to keywords.
% 
%           If the names of required and optional parameters were given:
%          Structure with fieldnames corresponding to the parameter names
%          and the values of those fields set to the parameter values.
%          Optional parameters that did not appear are set to their default 
%          values as given in the input argument par_opt.
%
% keyval    Structure with fieldnames corresponding to the keywords and
%          values that are read form the argument list, or from the default
%          values in keyval_def for those keywords that were not in the
%          argument list.
%
% present   Structure with field names matching the positional parameter names
%          (if they were given) and the keyword names, and which have values
%          logical 0 or 1 indicating if the parameter or keyword appeared in args.
%          If a keyword appeared as its negation e.g. 'nofoo', then it is deemed
%          to have been present i.e. present.foo = 1
%
% filled    Structure with field names matching the positional parameter names
%          (if they were given) and the keyword names, and which have values
%          logical 0 or 1 indicating if the argument is non-empty (whether
%          that be because it was supplied with a non-empty default, or
%          because it was given a non-empty value on the command line).
%
% ok        True if all is OK, false if not. If there is an error, but ok is
%          not a return argument, then an exception will be thrown. If ok is
%          a return argument, then an error will not throw an exception, so 
%          you must test the value of ok on return.
%
% mess      Error message of not ok; empty string if all is ok.
%
%
% EXAMPLE 1: Unlimited number of positional parameters:
% =====================================================
% Consider the function:
%       function parse_test (varargin)
%
%       % Argument names and default values:
%       keyval_def = struct('background',[12000,18000], ...
%                           'normalise', 1, ...
%                           'modulation', 0, ...
%                           'output', 'data.txt');
%
%       % Arguments which are logical flags:
%       flagnames = {'normalise','modulation'};
%
%       % Parse input:
%       [par, keyval, present] = symbz.parse_arguments...
%                   (varargin, keyval_def, flagnames);
%
%       % Display results
%       par
%       keyval
%       present
%
%       end
%
% Then calling parse_test with input as follows:
%   >> parse_test('input_file.dat',18,{'hello','tiger'},...
%                       'back',[15000,19000],'mod','nonorm')
%
% results in the output:
%   par =
%        'input_file.dat'    [18]    {1x2 cell}
%
%   argout =
%        background: [15000 19000]
%         normalise: 0
%        modulation: 1
%            output: 'data.txt'
%
%   present =
%        background: 1
%         normalise: 1
%        modulation: 1
%            output: 0
%
%
% EXAMPLE 2: Named required and optional parameters
% =================================================
%
%       function parse_test (varargin)
%
%       % Required parameters:
%       par_req = {'data_source', 'ei'};
%
%       % Optional parameters:
%       par_opt = struct('emin', -0.3, 'emax', 0.95, 'de', 0.005);
%
%       % Argument names and default values:
%       keyval_def = struct('background',[12000,18000], ...
%                           'normalise', 1, ...
%                           'modulation', 0, ...
%                           'output', 'data.txt');
%
%       % Arguments which are logical flags:
%       flagnames = {'normalise','modulation'};
%
%       % Parse input:
%       [par,argout,present] = symbz.parse_arguments(varargin,...
%                       par_req, par_opt, keyval_def, flagnames);
%
%       % Display results
%       par
%       argout
%       present
%
%       end
%
% Then calling parse_test with input:
%   >> parse_test('input_file.dat',18,-0.5,0.6,...
%                       'back',[15000,19000],'mod','nonorm')
%
% results in the output:
%     par = 
% 
%         data_source: 'input_file.dat'
%                  ei: 18
%                emin: -0.5000
%                emax: 0.6000
%                  de: 0.0050
% 
% 
%     argout = 
% 
%         background: [15000 19000]
%          normalise: 0
%         modulation: 1
%             output: 'data.txt'
% 
% 
%     present = 
% 
%         data_source: 1
%                  ei: 1
%                emin: 1
%                emax: 1
%                  de: 0
%          background: 1
%           normalise: 1
%          modulation: 1
%              output: 0


% Original author: T.G.Perring 
% 
% $Revision:: 830 ($Date:: 2019-04-09 10:03:50 +0100 (Tue, 9 Apr 2019) $) 
%

throw_error=(nargout<=4);

% Check if legacy input arguments format
if nargin>=2 && nargin<=4 && isstruct(varargin{1})
    % Legacy input: (args, keyval_def, [flagnames], [opt] )
    par_req=0;
    par_opt_def=Inf;
    keyval_def=varargin{1};
    nopt=2;
elseif nargin>=4
    % New format: (args,par_req,par_opt_def,keyval_def,varargin)
    par_req=varargin{1};
    par_opt_def=varargin{2};
    keyval_def=varargin{3};
    nopt=4;
else
    % Error:
    ok=false;
    mess='Check number and type of input arguments';
    if throw_error, error(mess), else...
            [par,keyval,present,filled]=error_return; return
    end
end

% Get flagnames and opt
[ok, mess, flagnames, opt] = parse_opt_args (varargin{nopt:end});
if ~ok
    if throw_error, error(mess), else...
            [par,keyval,present,filled]=error_return; return
    end
end

% Check parameter arguments
[ok,mess,par,nam_par,np_req,np_opt]=check_par_arguments(par_req,par_opt_def);
if ~ok
    if throw_error, error(mess), else...
            [par,keyval,present,filled]=error_return; return
    end
end

% Check keyword arguments and get variables for parsing
[ok,mess,keyval,nam,namchk,ind,flag,negflag]=check_keywords...
    (keyval_def,flagnames,opt);
if ~ok
    if throw_error, error(mess), else...
            [par,keyval,present,filled]=error_return; return
    end
end

% Check parameter names and keywords are collectively unique
if ~(isempty(nam_par) || isempty(nam))
    tmp=sort([nam_par;nam]);
    for i=2:numel(tmp)
        if strcmpi(tmp{i-1},tmp{i})
            ok=false;
            mess=['The combined parameter and keyword name list contains ''',...
                tmp{i},''' at least twice (with keyword negation &/or prefix.'];
            if throw_error, error(mess), else...
                    [par,keyval,present,filled]=error_return; return
            end
        end
    end
end

% Parse the input argument list
[ok,mess,par,keyval,present]=parse_args(args,par,nam_par,np_req,np_opt,...
    keyval,nam,namchk,ind,flag,negflag,opt);

if ok
    filled=false(numel(nam),1);
    for i=1:numel(filled)
        if ~isempty(keyval.(nam{i}))
            filled(i)=true;
        end
    end
    filled=cell2struct(num2cell(filled),nam);
else
    if throw_error, error(mess), else...
            [par,keyval,present,filled]=error_return; return
    end
end


%----------------------------------------------------------------------------------------
function [ok, mess, flagnames, opt] = parse_opt_args (varargin)
% Get optional arguments

narg=numel(varargin);
if narg==2 && (iscellstr(varargin{1}) || isempty(varargin{1})) && isstruct(varargin{2})
    flagnames = varargin{1};
    [ok,mess,opt] = update_opt (varargin{2});
    
elseif narg==1 && iscellstr(varargin{1})
    flagnames=varargin{1};
    [ok,mess,opt] = update_opt (struct());      % loads default
    
elseif narg==1 && isstruct(varargin{1})
    flagnames={};
    [ok,mess,opt] = update_opt (varargin{1});
    
elseif narg==0
    flagnames={};
    [ok,mess,opt] = update_opt(struct());   % loads default
    
else
    ok=false;
    mess='Check validity of optional arguments ''flagnames'' and/or ''opt''';
    flagnames={};
    opt=struct();
end


%----------------------------------------------------------------------------------------
function [ok,mess,opt_out] = update_opt (opt)
% Update the values of the fields in the options structure, checking also
% that the input structure does not contain unexpected fields (which will
% likely be due to typos, given that this is for developers, not users)
%
% Returns default option values if there is an error

ok=true;
mess='';

nam={'prefix';'prefix_req';'flags_noneg';'flags_noval';...
    'keys_exact';'keys_at_end';'keys_once';'noffset'};
val={'';true;false;false;false;true;true;0};

optnam=fieldnames(opt);
for i=1:numel(optnam)
    ind=find(strcmpi(optnam{i},nam),1);
    if ~isempty(ind)
        if ind==1 && ischar(opt.(optnam{i}))
            val{1}=opt.(optnam{i});     % this way we dont worry about letter case
            
        elseif ind==2 && islognumscalar(opt.(optnam{i}))
            val{2}=logical(opt.(optnam{i}));
            
        elseif ind==3 && islognumscalar(opt.(optnam{i}))
            val{3}=logical(opt.(optnam{i}));
            
        elseif ind==4 && islognumscalar(opt.(optnam{i}))
            val{4}=logical(opt.(optnam{i}));
            
        elseif ind==5 && islognumscalar(opt.(optnam{i}))
            val{5}=logical(opt.(optnam{i}));
            
        elseif ind==6 && islognumscalar(opt.(optnam{i}))
            val{6}=logical(opt.(optnam{i}));
            
        elseif ind==7 && islognumscalar(opt.(optnam{i}))
            val{7}=logical(opt.(optnam{i}));
            
        elseif ind==8 && isnumeric(opt.(optnam{i})) && opt.(optnam{i})>=0
            val{8}=logical(opt.(optnam{i}));
        end
        
    else
        ok=false;
        mess=['Unexpected option name, type or value: ''',optnam{i},''''];
        
    end
end

opt_out=cell2struct(val,nam);


%----------------------------------------------------------------------------------------
function [ok,mess,par,nam,nreq,nopt]=check_par_arguments(par_req,par_opt)
% Check the validity of the arguments defining required and positional parameters
%
% Input:
% ------
%   par_req Required parameters: a cell array of parameter names or a number
%           0,1,2,...
%
%   par_opt Optional parameters: a structure giving names and default values
%           or a number 0,1,2,... or Inf
%
% In either case, an empty argument means no parameters of that type are
% permitted i.e. is equivalent to 0.
%
% Output:
% -------
%   ok      True if all OK, false otherwise
%
%   par     Structure with default values if parameters are named. (Note:
%          required parameters are given the value [], but by definition
%          these will be required to be given a value.)
%           Empty cell array if parameters are not named.
%
%   nam     Names of parameters. Empty cell array if parameters are not named.
%
%   nreq    Number of required parameters: 0,1,2,...
%
%   nopt    Number of optional parameters: 0,1,2,...  or Inf



% Check required parameter names
%-------------------------------
if ~isempty(par_req)
    if isnumeric(par_req)
        if isscalar(par_req) && par_req>=0 && isfinite(par_req) && rem(par_req,1)==0
            nreq=par_req;
            nam_req={};
        else
            ok=false;
            mess='The number of required parameters must be 0,1,2,...';
            par={}; nam={}; nreq=0; nopt=0;
            return
        end
        
    elseif iscellstr(par_req)
        for i=1:numel(par_req)
            if ~isvarname(par_req{i})
                ok=false;
                mess=['Required parameter name ''',par_req{i},...
                    ''' is not a valid Matlab variable name'];
                par={}; nam={}; nreq=0; nopt=0;
                return
            end
        end
        if numel(par_req)>1
            tmp=sort(par_req);
            for i=2:numel(par_req)
                if strcmpi(tmp{i-1},tmp{i})
                    ok=false;
                    mess=['The required parameters name list contains ''',...
                        tmp{i},''' at least twice. Names must be unique.'];
                    par={}; nam={}; nreq=0; nopt=0;
                    return
                end
            end
        end
        nreq=numel(par_req);
        nam_req=par_req(:);     % ensure column vector
    else
        ok=false;
        mess=['Required parameters must be given as a list of names or ',...
            'the number of required parameters'];
        par={}; nam={}; nreq=0; nopt=0;
        return
    end
else
    % Empty argument assumed to mean no required parameters
    nreq=0;
    nam_req={};
end


% Optional parameters
%--------------------
if ~isempty(par_opt)
    if isnumeric(par_opt)
        if isscalar(par_opt) && par_opt>=0 && (isinf(par_opt) || rem(par_opt,1)==0)
            nopt=par_opt;
            nam_opt={};
        else
            ok=false;
            mess=['The number of optional parameters must be 0,1,2,...',...
                ', or Inf for an unlimited number'];
            par={}; nam={}; nreq=0; nopt=0;
            return
        end
    elseif isstruct(par_opt)
        nam_opt=fieldnames(par_opt);
        if numel(nam_opt)>0
            nopt=numel(nam_opt);
        else
            nopt=0;
            nam_opt={};
        end
    else
        ok=false;
        mess=['Optional parameters must given as a structure or the ',...
            'number of required parameters'];
        par={}; nam={}; nreq=0; nopt=0;
        return
    end
else
    % Empty argument assumed to mean an unlimited number of optional parameters
    nopt=0;
    nam_opt={};
end


% Combine the input
% -----------------
% If required and optional parameters are both >0, must both be un-named or
% both named.
if nreq>0 && nopt>0
    % Both required parameters and optional parameters
    if ~xor(isempty(nam_req),isempty(nam_opt))
        if ~isempty(nam_req)     % both defined by names
            nam=[nam_req;nam_opt];
            tmp=sort(nam);
            for i=2:numel(tmp)
                if strcmpi(tmp{i-1},tmp{i})
                    ok=false;
                    mess=['The combined required and optional parameter names contain ''',...
                        tmp{i},''' at least twice. Names must be unique.'];
                    par={}; nam={}; nreq=0; nopt=0;
                    return
                end
            end
            par=catstruct(cell2struct(repmat({[]},nreq,1),nam_req), par_opt);
        else
            par={};
            nam={};
        end
    else
        ok=false;
        mess=['If a non-zero number of each of required and optional parameters',...
            ' is given, they must both be given numerically or both named'];
        par={}; nam={}; nreq=0; nopt=0;
        return
    end
elseif nreq>0
    % Required parameters only
    if ~isempty(nam_req)
        par=cell2struct(repmat({[]},nreq,1),nam_req);
        nam=nam_req;
    else
        par={};
        nam={};
    end
elseif nopt>0
    % Optional parameters only
    if ~isempty(nam_opt)
        par=par_opt;
        nam=nam_opt;
    else
        par={};
        nam={};
    end
else
    % No required and no optional parameters
    par={};
    nam={};
end


% Ok if got to here
%------------------
ok=true;
mess='';


%----------------------------------------------------------------------------------------
function [ok,mess,keyval,nam,namchk,ind,flag,negflag]=check_keywords...
    (keyval_in,flagnam,opt)
%
%
% Input:
% ------
%   keyval_in   Structure: fieldnames are names of keywords, field values
%              are the default values
%   flagnam     Call array with names of flags. The names must fields of
%              keyval_in
%   opt         Structure with values of options: required fields are:
%                   prefix      Value of prefix to keywords. Must be a
%                              character string.
%                   prefix_req  True if require that the prefix be present,
%                              false otherwise
%                   flags_noneg True if flags cannot be negated.
%                   flags_noval True if flags cannot be given values. If
%                              flags_noneg and flags_noval then the default
%                              for flags must be false, as there is no way
%                              to set them to false.
%
% Output:
% -------
%   ok          true if no problems constructing keyword list
%   mess        '' id ok, otherwise contains error message
%   keyval      Structure with defaults for flags turned in logicals
%   nam         Keywords i.e. fieldnames(keyval)
%   namchk      Keywords with all permissible prefixing or negation. This will
%              be used to parse arguments. Each keyword is unique.
%   ind         Index of permissible keyworwds into the root keyword list
%   flag        Logical array with true where namchk is a flag or negation
%              of a flag
%   negflag     Logical array with true where namchk is the negation of a flag
%
%
% If a keyword='hello' and prefix='-', then if opt.prefix_req==false valid
% keywords are:
%       'hello' and '-hello'
%
% If the keyword is also a flag, then the full list is:
%       'hello', '-hello', 'nohello' and '-nohello'
%
% Checka re performed to ensure that prefix and flag status do not construct
% duplicate names. For example:
% - 'nohello' and 'hello' are not permitted if 'hello' is a flag
%       (because the negation of 'hello' matches 'nohello').
% - 'nohello' and 'hello' are permitted if 'hello' is NOT a flag


% Note about efficiency: TGP did some timing tests of the 'N^2' algorithm
% where for each keyword we check against all others, and an algorithm where
% the keywords are sorted and then compare adjacent keywords. For 300 keywords
% or less, the N^2 algorithm is about 5 times slower, and about 2x slower
% for N<20. Given the fact that in practice N<20, the other overheads,
% that the calling function will march through the list anyway at least once,
% it is not worth the effort of optimising where the comparison keyword comes
% from a second list that also needs to be incremented.


% Catch case of no keywords
if isempty(keyval_in)
    % Any empty object
    keyval=struct([]);
    nam=cell(0,1);
    
elseif isstruct(keyval_in)
    % Keyword structure given
    keyval=keyval_in;
    nam=fieldnames(keyval);
    
else
    keyval=struct([]);
    nam=cell(0,1);
    n=0;
    flag=false(0,1);
    ok=false;
    mess='Keywords and their defaults can only be given as a structure';
    namchk=nam; ind=(1:n)'; negflag=false(n,1);
    return
end

n=numel(nam);
flag=false(n,1);

% Determine which names are flags
if ~isempty(flagnam)
    for i=1:numel(flagnam)
        ipos = find(strcmpi(flagnam{i},nam),1);    % can only be 0 or 1 match
        if ~isempty(ipos)
            flag(ipos)=true;
            % Ensure default is a logical scalar
            val=keyval.(nam{ipos});
            if islognumscalar(val)
                if ~logical(val) || ~opt.flags_noneg || ~opt.flags_noval
                    keyval.(nam{ipos}) = logical(val);
                else
                    ok=false;
                    mess=['Default value of flag ''',nam{ipos},...
                        ''' must be false if flags_noneg and flags_noval are both true'];
                    namchk=nam; ind=(1:n)'; negflag=false(n,1);
                    return
                end
            else
                ok=false;
                mess=['Default value of flag ''',nam{ipos},...
                    ''' must be 0 or 1, or true or false'];
                namchk=nam; ind=(1:n)'; negflag=false(n,1);
                return
            end
        else
            ok=false;
            mess=['Flag name ''',flagnam{i},''' is not in the named parameter list'];
            namchk=nam; ind=(1:n)'; negflag=false(n,1);
            return
        end
    end
end

% Determine full list of names allowing for prefix and negation of flags
if any(flag) && ~opt.flags_noneg
    negnam=nam(flag);
    for i=1:numel(negnam)
        negnam{i}=['no',negnam{i}];
    end
    namchk=[nam;negnam];
    ind=[(1:n)';find(flag)];
    flag=[flag;true(numel(negnam),1)];
    negflag=[false(n,1);true(numel(negnam),1)];
else
    namchk=nam;
    ind=(1:n)';
    negflag=false(n,1);
end

if ~isempty(opt.prefix)
    if opt.prefix_req
        for i=1:numel(namchk)
            namchk{i}=[opt.prefix,namchk{i}];
        end
    else
        prenamchk=cell(size(namchk));
        for i=1:numel(namchk)
            prenamchk{i}=[opt.prefix,namchk{i}];
        end
        namchk=[namchk;prenamchk];
        ind=[ind;ind];
        flag=[flag;flag];
        negflag=[negflag;negflag];
    end
end

% Check there are no duplications in the list arising from negation and/or prefixing
if numel(namchk)>1
    [tmp,ix]=sort(namchk);
    for i=2:numel(namchk)
        if strcmpi(tmp{i-1},tmp{i})
            ok=false;
            mess=['The keywords ''',nam{ind(ix(i-1))},''' and ''',nam{ind(ix(i))},...
                ''' are ambiguous due to prefixing and/or negation'];
            return
        end
    end
end

% All ok if got to here
ok=true;
mess='';


%----------------------------------------------------------------------------------------
function [ok,mess,par,keyval,present]=parse_args(args,par,nam_par,np_req,np_opt,keyval,...
    nam,namchk,ind,flag,negflag,opt)
% Parse the input argument list
%
%   opt     Options structure required fields
%               flags_noval
%               keys_exact
%               keys_at_end
%               keys_once
%               noffset

ok=true;
mess='';

narg=numel(args);
indpar=false(1,narg);
np_max=np_req+np_opt;
nkey=numel(nam);

i=1;
np=0;
key_present=false(numel(nam),1);
expect_key=false;
while i<=narg
    % Determine if argument is a keyword; ambiguous keywords are an error
    iskey=false;
    if nkey>0 && ischar(args{i}) && ~isempty(args{i})
        if opt.keys_exact
            ipos=find(strcmpi(args{i},namchk)); % strings in namchk are unique
            if ~isempty(ipos), iskey=true; end
        else
            ipos=stringmatchi(args{i},namchk);
            if ~isempty(ipos)
                if numel(ipos)==1
                    iskey=true;
                else
                    ok=false;
                    mess=['Ambiguous keyword at position ',...
                        num2str(i+opt.noffset),' in the argument list'];
                    break
                end
            end
        end
    end
    
    % Branch on parameter or keyword
    if ~iskey
        if ~expect_key
            np=np+1;
            if np<=np_max
                indpar(i)=true;
            else
                ok=false;
                mess=['The number of positional parameter(s) exceeds the maximum request of ',num2str(np_max)];
                break
            end
        else
            ok=false;
            mess=['Expected a keyword but found ',disp_string(args{i}),...
                ' at position ',num2str(i+opt.noffset),' in the argument list'];
            break
        end
        
    else
        ikey=ind(ipos);
        % Check if the keyword has already appeared or not
        if ~key_present(ikey)
            key_present(ikey)=true;
        elseif opt.keys_once
            ok=false;
            mess=['Keyword ''',nam{ikey},...
                ''' (or its negation if a flag) appears more than once'];
            break
        end
        % Get value corresponding to keyword
        if flag(ipos)
            if opt.flags_noval || i==narg || ~islognumscalar(args{i+1})
                keyval.(nam{ikey})=~negflag(ipos);
            else
                i=i+1;
                keyval.(nam{ikey})=xor(negflag(ipos),logical(args{i}));
            end
        else
            if i<narg
                i=i+1;
                keyval.(nam{ikey})=args{i};
            else
                ok=false;
                mess=['Keyword ''',nam{ikey},...
                    ''' expects a value, but the keyword is the final argument'];
                break
            end
        end
        expect_key = opt.keys_at_end;   % if true, will expect only keywords from now
    end
    i=i+1;
end

% Searched the argument list. Now pack final output
if ok
    if np>=np_req
        if isempty(par) % no named parameters
            par=args(indpar);
            present=cell2struct(num2cell(key_present),nam);
        else
            if np>0
                ix=find(indpar);
                for i=1:np
                    par.(nam_par{i})=args{ix(i)};
                end
                par_present=[true(np,1);false(np_max-np,1)];
                present=cell2struct(num2cell([par_present;key_present]),[nam_par;nam]);
            else
                present=cell2struct(num2cell(key_present),nam);
            end
        end
    else
        ok=false;
        mess=['The number of positional parameter(s) is less than the minimum request of ',num2str(np_req)];
        present=struct([]);
    end
else
    present=struct([]);
end


%----------------------------------------------------------------------------------------
function str=disp_string(var)
% Encapsulate information about enexpected variable as a string

nchar_max=10;
if ischar(var)
    if isempty(var)
        str='empty string';
    else
        if numel(var)<nchar_max
            str=['string ''',var,''''];
        else
            str=['string ''',var(1,nchar_max),'...'''];
        end
    end
else
    if isscalar(var)
        str=['argument of class ',class(var)];
    else
        sz=size(var);
        sz_str='[';
        for i=1:numel(sz)
            sz_str=[sz_str,num2str(sz(i)),'x'];
        end
        sz_str(end)=']';
        str=['argument ',sz_str,' array of class ',class(var)];
    end
end


%----------------------------------------------------------------------------------------
function [par,keyval,present,filled]=error_return
% Default output if there is an error
par={};
keyval=struct([]);
present=struct([]);
filled=struct([]);

function ok=islognumscalar(val)
% Determine if a value is a (non-empty) scalar logical, or is numeric 0 or 1
if isscalar(val) && (islogical(val) || (isnumeric(val) && (val==0 ||val==1)))
    ok=true;
else
    ok=false;
end

function ind = stringmatchi (str, strcell)
% Index of match(es) or unambiguous abbreviations in a cell array of strings
%
%   >> ind = stringmatchi (string,strcell)
%
% Input:
% ------
%   string  Test string
%   strcell Cell array of strings
%
% Output:
% -------
%   ind     Index of str in strcell if str is an exact match or unambiguous
%          abbreviation of one of the elements of cellst.
%           If str is an exact match for one or more elements of strcell,
%          only these indicies are returned even if it is also an abbreviation
%          of other element


if ~ischar(str)
    error('First argument must be a string')
end

nch=numel(str);
ind=find(strncmpi(str,strcell,nch));

% If string and cellstr and more than one match, look for equality
if numel(ind)>1
    ix=false(size(ind));
    for i=1:numel(ind(:))
        if numel(strcell{ind(i)})==nch
            ix(i)=true;
        end
    end
    if any(ix(:))
        ind=ind(ix);
    end
end