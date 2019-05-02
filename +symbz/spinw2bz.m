function bz = spinw2bz(sw,varargin)
[~,rlat]=symbz.spinw2lat(sw,varargin);
bz = symbz.brillouinzone(rlat,varargin);
end