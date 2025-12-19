function [rch,gch,bch,caO]=srhist_color_ca(sz,zm,xtot,ytot,ttot,cmp,ca)
if nargin<6
    cmp = colormap(hsv(64));
    ca = [];
elseif nargin<7
    ca = [];
end
ttot=ttot(1:size(xtot,1));

segnum = size(cmp,1);

if isempty(ca)
    Mt = max(ttot);
    mt = min(ttot);
else
    Mt = ca(2);
    mt = ca(1);
end

caO = [mt,Mt];

incre=floor((Mt-mt+1)/segnum);
for ii=1:1:segnum
    tst=(ii-1)*incre + mt;
    ted=ii*incre + mt;
    
    mask=ttot>=tst & ttot<=ted;
    xtmp=xtot(mask);
    ytmp=ytot(mask);
    
    tmpim=SRreconstructhist(sz,zm,xtmp,ytmp);
    tmpim=double(tmpim);
    if ii==1
        rch=cmp(ii,1)*tmpim;
        gch=cmp(ii,2)*tmpim;
        bch=cmp(ii,3)*tmpim;
    else
        rch=rch+cmp(ii,1)*tmpim;
        gch=gch+cmp(ii,2)*tmpim;
        bch=bch+cmp(ii,3)*tmpim;
    end
end
