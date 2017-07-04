% fit ZXY cell-array scatter points to ellipsoid of minimum rms (affine transformed) unit-sphere width
%
% INPUT
%   zxy: Nx1 cell array of ZXY data points
%   param0: [[CENT_XYZ_0],[EUL_ANGLE_0],[ERAD_0]]
% OUTPUT
%   param_ellipsoid: best ellipsoid param found
%   w_usph: rms width of unit-sph affine transformed from ellipsoid params

function [param_ellipsoid,w_usph] = fit_ellipsoid(zxy,param0)
    halo_zxy=vertcat(zxy{:});   % concatenate all shots into single array
    
    % get initial param guess
    cent0=param0(1:3);
    eeul0=param0(4:6);
    erad0=param0(7:9);
    
    % build parameter constraints
    centlb=cent0-erad0;
    eeullb=-pi*[1,1/2,1];
    eradlb=erad0/2;
    paramlb=[centlb,eeullb,eradlb];
    
    centub=cent0+erad0;
    eeulub=pi*[1,1/2,1];
    eradub=erad0*2;
    paramub=[centub,eeulub,eradub];
    
    % define optimisation problem
    tol=1e-12;
    foption=optimset('Display','iter-detailed','FunValCheck','on',...
        'MaxFunEvals',1e5,'MaxIter',1e5,...
        'TolFun',tol,'TolX',tol,...
        'PlotFcns',{@optimplotx,@optimplotfval,@optimplotfunccount});
    foption.StepTolerance=tol;
    foption.FunctionTolerance=tol;
    
    A=[];
    b=[];
    Aeq=[];
    beq=[];
    nonlcon=[];
    
    % run the optimiser!
    [param_ellipsoid,w_usph,exitflag,output]=fmincon(@(x)halo_sph_rms(halo_zxy,x),...
        param0,A,b,Aeq,beq,paramlb,paramub,nonlcon,foption);
end