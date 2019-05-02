function bzg = spinw2bzg(sw,varargin)
[~,rlat]=symbz.spinw2lat(sw,varargin);
bz = symbz.brillouinzone(rlat,varargin);
bzg = symbz.BZGridQ(bz,varargin{:},'complex',true);
end