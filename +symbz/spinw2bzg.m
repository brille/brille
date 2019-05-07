function [bzg,trnm] = spinw2bzg(sw,varargin)
[~,rlat,trnm]=symbz.spinw2lat(sw,varargin);
bz = symbz.brillouinzone(rlat,varargin);
bzg = symbz.BZGridQ(bz,varargin{:},'complex',true);
end