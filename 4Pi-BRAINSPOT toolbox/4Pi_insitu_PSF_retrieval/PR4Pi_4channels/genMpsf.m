%% Generate normalized measured PSF that are used for phase retrieval.

%%
function obj = genMpsf(obj)

R1=obj.SubroiSize;    
R=obj.PSFsize;
N=obj.DatadimZ;
N_ch = obj.DatadimCh;

[X1,Y1]=meshgrid(-R1/2:R1/2-1,-R1/2:R1/2-1);
circle_tmp=sqrt(X1.^2+Y1.^2);
circleF=circle_tmp;
circleF(circleF<=R1/2-1)=1;
circleF(circleF>R1/2-1)=0;
if R==R1
    circleF=ones(R,R);
end
realsize0=floor(R1/2);
realsize1=ceil(R1/2);
starty=-realsize0+R/2+1;endy=realsize1+R/2;
startx=-realsize0+R/2+1;endx=realsize1+R/2;
MpsfL=zeros(R,R,N,N_ch); 
% mask index
index_mask1 = circle_tmp<=R1/2-1 & circle_tmp>=R1/2-3 & X1 >0 & Y1 > 0;
index_mask2 = circle_tmp<=R1/2-1 & circle_tmp>=R1/2-3 & X1 <=0 & Y1 > 0;
index_mask3 = circle_tmp<=R1/2-1 & circle_tmp>=R1/2-3 & X1 >0 & Y1 <= 0;
index_mask4 = circle_tmp<=R1/2-1 & circle_tmp>=R1/2-3 & X1 <=0 & Y1 <= 0;
for jj = 1 : N_ch 
    for ii=1:N
        Mpsfo=obj.Mpsf_subroi(:,:,ii,jj);  

        Edge=[mean(Mpsfo(index_mask1))+1.*std(Mpsfo(index_mask1)),mean(Mpsfo(index_mask2))+1.*std(Mpsfo(index_mask2)),...
            mean(Mpsfo(index_mask3))+1.*std(Mpsfo(index_mask3)),mean(Mpsfo(index_mask4))+1.*std(Mpsfo(index_mask4))];  
 
        bg=max(Edge);

        Fig2=(Mpsfo-bg);
        Fig2 = Fig2 .* circleF;
        Fig2(Fig2<=0)=0;
        tmp=zeros(R,R);
        tmp(starty:endy,startx:endx)=Fig2;
        Fig2=tmp;
        minimum0=0;
        Fig2(Fig2<=minimum0)=minimum0;
        MpsfL(:,:,ii,jj)=Fig2./sum(sum(Fig2));
        
        Fig3=Mpsfo;
        Fig3 = Fig3 .* circleF;
        Fig3(Fig3<=0)=0;
        tmp=zeros(R,R);
        tmp(starty:endy,startx:endx)=Fig3;
        Fig3=tmp;
        minimum0=0;
        Fig3(Fig3<=minimum0)=minimum0;
        MpsfL_2(:,:,ii,jj)=Fig3./sum(sum(Fig3));
    end
end

obj.Mpsf_extend = MpsfL;
obj.Mpsf_extend_v2 = MpsfL_2;





