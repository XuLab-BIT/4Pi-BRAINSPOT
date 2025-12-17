%% precomputeParam - generate images for k space operation, and saved in precomputed parameters.

%%            
function obj = precomputeParam(obj)

obj.Mpsf_subroi = obj.BeadData;
[obj.DatadimY,obj.DatadimX,obj.DatadimZ,obj.DatadimCh]=size(obj.Mpsf_subroi);

[XC,YC]=meshgrid(-obj.DatadimX/2:obj.DatadimX/2-1,-obj.DatadimY/2:obj.DatadimY/2-1);
obj.PhiC=atan2(YC,XC);
obj.ZoC=sqrt(XC.^2+YC.^2);

[X,Y]=meshgrid(-obj.PSFsize/2:obj.PSFsize/2-1,-obj.PSFsize/2:obj.PSFsize/2-1);
obj.Zo=sqrt(X.^2+Y.^2);
scale=obj.PSFsize*obj.Pixelsize;
obj.k_r=obj.Zo./scale;
obj.Phi=atan2(Y,X);
Freq_max=obj.PRstruct.NA/obj.PRstruct.Lambda;
obj.NA_constrain=obj.k_r<Freq_max;

% create Zernike_Polynomials object
zk = Zernike_Polynomials();
zk.Ordering = 'Wyant';
zk.setN(obj.ZernikeorderN);
zk.initialize();
[Zrho, Ztheta, Zinit] = ...
    zk.params3_Zernike(obj.Phi, obj.k_r, obj.PRstruct.NA, obj.PRstruct.Lambda);
zk.matrix_Z(Zrho, Ztheta, Zinit);
obj.Z = zk;

