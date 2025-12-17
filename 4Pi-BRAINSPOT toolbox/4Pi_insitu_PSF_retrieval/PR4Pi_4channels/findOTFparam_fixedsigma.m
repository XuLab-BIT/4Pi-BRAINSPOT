function obj=findOTFparam_fixedsigma(obj,sigma)
% findOTFparam - fixed SigmaX and SigmaY of a Gaussian filter for OTF rescale. 
R=obj.PSFsize;
z = obj.Zpos;
N=obj.DatadimZ;
N_ch = obj.DatadimCh;


% Calculate the otf scale using summed PSFs
if nargin == 1
    [val,ind]=min(abs(z));
    
    realsize0=floor(obj.OTFratioSize/2);
    realsize1=ceil(obj.OTFratioSize/2);
    starty=-realsize0+R/2+1;endy=realsize1+R/2;
    startx=-realsize0+R/2+1;endx=realsize1+R/2;

    Mpsf_sum = squeeze(sum(sum(obj.Mpsf_extend,3),4));
    PRpsf_sum = squeeze(sum(sum(obj.PSFstruct.PRpsf,3),4));
    
    
    mOTF=fftshift(ifft2(Mpsf_sum));% measured OTF
    rOTF=fftshift(ifft2(PRpsf_sum)); % phase retrieved OTF
    tmp=abs(mOTF)./abs(rOTF);
    tmp1=tmp(startx:endx,starty:endy);
    ratio=tmp1;

    fit_param=[2,2,0];
    [I,sigma,bg,fit_im]=GaussRfit(obj,fit_param,ratio);
    


else  %using default sigma
    
    scale=R*obj.Pixelsize;
    [xx,yy]=meshgrid(-R/2:R/2-1,-R/2:R/2-1);
    X=abs(xx)./scale;
    Y=abs(yy)./scale;
    fit_im=1.*exp(-X.^2./2./sigma^2).*exp(-Y.^2./2./sigma^2);
end


% generate zernike fitted PSF modified by OTF rescale
Mod_psf=zeros(R,R,N,N_ch);
for jj = 1:N_ch
    for ii=1:N
        Fig4=obj.PSFstruct.PRpsf(:,:,ii,jj);   
        Fig4=Fig4./sum(sum(Fig4));
        Mod_OTF=fftshift(ifft2(Fig4)).*fit_im;
        Fig5=abs(fft2(Mod_OTF));
        Mod_psf(:,:,ii,jj)=Fig5;
    end
end
obj.PRstruct.SigmaX=sigma;
obj.PRstruct.SigmaY=sigma;
obj.PSFstruct.Modpsf=Mod_psf;
end


function [I,sigmax,sigmay,bg]=GaussRfit_v2(obj,ratio)
R=obj.PSFsize;
R1=obj.OTFratioSize;
scale=R*obj.Pixelsize;
x = [-R1/2:R1/2-1]./scale;
Ix = mean(ratio,1);
Iy = mean(ratio,2);

fx = fit(x',Ix','gauss1');
sigmax = fx.c1/sqrt(2);
fy = fit(x',Iy,'gauss1');
sigmay = fy.c1/sqrt(2);

I = (fx.a1 + fy.a1)/2;

Ixf = feval(fx,x);
Iyf = feval(fy,x);
figure;plot(x,Ix,'bo',x,Ixf,'r-')
figure;plot(x,Iy,'bo',x,Iyf,'r-')
bg = 0;

end

function [fit_im] = genfitim(obj,OTFparams)
R=obj.PSFsize;
I = OTFparams(1);
sigmax = OTFparams(2);
sigmay = OTFparams(3);
bg=0;
scale=R*obj.Pixelsize;
[xx,yy]=meshgrid(-R/2:R/2-1,-R/2:R/2-1);

X=abs(xx)./scale;
Y=abs(yy)./scale;
fit_im=I.*exp(-X.^2./2./sigmax^2).*exp(-Y.^2./2./sigmay^2)+bg;
end

function [I,sigma,bg,fit_im]=GaussRfit(obj,startpoint,input_im)

R1=obj.OTFratioSize;
R=obj.PSFsize;
estimate=fminsearch(@(x) Gauss2(x,input_im,R1,R,obj.Pixelsize),startpoint,optimset('MaxIter',50,'Display','off'));

I=estimate(:,1);
tmp=estimate(:,2);
tmp(tmp>5)=5;
sigma = tmp;

bg=0;

scale=R*obj.Pixelsize;
[xx,yy]=meshgrid(-R/2:R/2-1,-R/2:R/2-1);

X=abs(xx)./scale;
Y=abs(yy)./scale;
fit_im=I.*exp(-X.^2./2./sigma^2).*exp(-Y.^2./2./sigma^2)+bg;
end


function [sse,Model]=Gauss2(x,input_im,R1,R,pixelsize)
I=x(1);

sigma=x(2);
bg=0;
[xx,yy]=meshgrid(-R1/2:R1/2-1,-R1/2:R1/2-1);

scale=R*pixelsize;
X=abs(xx)./scale;
Y=abs(yy)./scale;
Model=I.*exp(-X.^2./2./sigma^2).*exp(-Y.^2./2./sigma^2)+bg;

sse=sum(sum((Model-input_im).^2));
end