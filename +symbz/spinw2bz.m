function [bz,trnm] = spinw2bz(sw,varargin)
[~,rlat,trnm]=symbz.spinw2lat(sw,varargin);
bz = symbz.brillouinzone(rlat,varargin);
end