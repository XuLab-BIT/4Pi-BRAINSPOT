function obj = datapreprocess(obj)

R1 = obj.SubroiSize;
N = obj.DatadimZ;
N_ch = obj.DatadimCh;

realsize0=floor(R1/2);
realsize1=ceil(R1/2);
MpsfC=zeros(R1,R1,N,N_ch);
shiftxy=obj.Beadcenter-obj.BeadXYshift;
for jj = 1 : N_ch 
    for ii=1:N
        tmp=fftshift(ifft2(squeeze(obj.BeadData(:,:,ii,jj))));
        shiftphase=-obj.ZoC./obj.DatadimX.*cos(obj.PhiC).*(obj.DatadimX/2-shiftxy(1)-1)-obj.ZoC./obj.DatadimY.*sin(obj.PhiC).*(obj.DatadimY/2-shiftxy(2)-1);
        tmp1=fft2(tmp.*exp(-2*pi.*shiftphase.*1i));
        startx=-realsize0+obj.DatadimX/2+1;endx=realsize1+obj.DatadimX/2;
        starty=-realsize0+obj.DatadimY/2+1;endy=realsize1+obj.DatadimY/2;
        tmp2=abs(tmp1(starty:endy,startx:endx));
        MpsfC(:,:,ii,jj)=tmp2;
    end
end
obj.Mpsf_subroi=MpsfC;
