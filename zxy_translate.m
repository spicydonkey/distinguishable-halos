function ZXY_OUT = zxy_translate(ZXY_IN,DISP,SIGN)
% Translates uniformly all ZXY-count of each shot by the correponding displacement vector
% DKS 1/11/16
%
% I/O
% ZXY_IN: should be a Nx1 cell of ZXY-array
% DISP: Nx1 cell of a single ZXY displacement vector
% SIGN: +/-1 adds/subtracts displacement (default=ADDITIVE)

if ~exist('SIGN','var')
    SIGN=1;     % default sign is ADDITION
end

nShot=size(ZXY_IN,1);

ZXY_OUT=cell(nShot,1);
for i=1:nShot
    nCount=size(ZXY_IN{i},1);
    ZXY_OUT{i}=ZXY_IN{i}+SIGN*repmat(DISP{i},[nCount,1]);
end