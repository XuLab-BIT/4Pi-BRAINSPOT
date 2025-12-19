function Analysis_spine_section_3D(obj,para,subname,tagS)

%%  Parameter

sz = para.sz;
psz = para.psz;
zm = para.zm;
sigma = para.sigma;
Recon_color_highb = para.Recon_color_highb*0.8;

tk = para.tk/psz;      %% convert to pixel    
wk = para.wk/psz;   %%  convert to pixel

cm = para.cm;
hh = sz/2+0.5;

savepath = para.savepath;

resultpath = [savepath, subname, '\'];
if ~exist(resultpath,'dir')
    mkdir(resultpath)
end
resultpath2 = [savepath, subname, '\Area\'];
if ~exist(resultpath2,'dir')
    mkdir(resultpath2)
end

reconx = obj.y./psz;
recony = obj.x./psz;
reconz = obj.z;

%%  Hand mark points for structural selection

if (tagS>=1 && exist(fullfile(resultpath,'Select.txt'),'file'))
    fileID = fopen(fullfile(resultpath,'Select.txt'),'r');
    formatSpec = '%d %d';
    sizeA = [2 Inf];
    xy_zm = fscanf(fileID,formatSpec,sizeA)';
    xy_zm = (xy_zm-1)/50*zm+1;
    fclose(fileID);

else
    ds = 4;
    [rch, gch, bch] = srhist_color_ca(sz,zm/ds,reconx,recony,reconz,cm);

    rchsm = gaussf(rch,sigma(1)/ds);
    gchsm = gaussf(gch,sigma(1)/ds);
    bchsm = gaussf(bch,sigma(1)/ds);

    rchsmst = imstretch_linear(rchsm,0,Recon_color_highb*ds,0,255);
    gchsmst = imstretch_linear(gchsm,0,Recon_color_highb*ds,0,255);
    bchsmst = imstretch_linear(bchsm,0,Recon_color_highb*ds,0,255);

    colorim = double(joinchannels('RGB',rchsmst,gchsmst,bchsmst));
    im = rchsmst+gchsmst+bchsmst;

    disp('Right click to stop picking points');
    figure('Position',[10,10,1000,1000])
    [xy_zm] = pickspot(im*1.5);
    xy_zm = xy_zm*ds+1;
end
close all
%%  Generate XY projection with selection curve

disp('Rending pseudo-colored imaging');
if isfield(para,'ca')
    [rch, gch, bch,ca0] = srhist_color_ca(sz,zm,reconx,recony,reconz,cm,para.ca);
else
    [rch, gch, bch,ca0] = srhist_color_ca(sz,zm,reconx,recony,reconz,cm);
end

    
rchsm = gaussf(rch,sigma(1));
gchsm = gaussf(gch,sigma(1));
bchsm = gaussf(bch,sigma(1));

rchsmst = imstretch_linear(rchsm,0,Recon_color_highb,0,255);
gchsmst = imstretch_linear(gchsm,0,Recon_color_highb,0,255);
bchsmst = imstretch_linear(bchsm,0,Recon_color_highb,0,255);

colorim = double(joinchannels('RGB',rchsmst,gchsmst,bchsmst));
im = rchsmst+gchsmst+bchsmst;

tmp = colorim(colorim>2);
[val,ed] = histcounts(tmp,100,'Normalization','cdf');
indct = find(val>0.99);
normf = ed(indct(1));

imwrite(colorim/normf,fullfile(resultpath,'Reconstruction_XY.tiff'));


ri = xy_zm(end,1)-xy_zm(1,1)+1i*xy_zm(end,2)-1i*xy_zm(1,2);
ra = angle(ri);

rTL = [cos(ra) sin(ra) 0; -sin(ra) cos(ra) 0; 0 0 1]*[xy_zm'-hh*zm;zeros(1,length(xy_zm))] +hh*zm;

fC = fit(rTL(1,:)',rTL(2,:)','smoothingspline');
rx_zm = linspace(rTL(1,1),rTL(1,end),10*length(rTL))';
ry_zm = feval(fC,rx_zm);
curveLength = sum(vecnorm(  diff( [rx_zm, ry_zm] )  ,2,2));
nn = floor(curveLength/tk/zm);

pT = [cos(ra) -sin(ra) 0; sin(ra) cos(ra) 0; 0 0 1]*[rx_zm'-hh*zm;ry_zm'-hh*zm;zeros(size(rx_zm'))] +hh*zm;
pt = interparc(nn,pT(1,:)',pT(2,:)');
XX_zm = pt(:,1);
YY_zm = pt(:,2);
XX = XX_zm/zm;
YY = YY_zm/zm;

xy_zm = (xy_zm-1)/zm*50+1;
fileID = fopen(fullfile(resultpath,'Select.txt'),'w');
fprintf(fileID,'%10d %10d\n',round(xy_zm'));
fclose(fileID);

%%  Generate cross sections along the curve

h = figure('visible','on','MenuBar','none','Toolbar','none','Name','Analysis ROI','Position',[0,0,1150,1150]);
hold on
h.InvertHardcopy = 'off';
ax = axes('Position',[0,0,1,1]);
image(colorim./normf);
axis equal 
axis off
line(XX_zm,YY_zm,'color','y','linew',0.75)
print(gcf,fullfile(resultpath,'Selected_ROI_1'),'-dpng','-r300');


App = struct;
roi = zeros(1,4);
for i = 1:length(XX)-1
    ri = XX(i+1)-XX(i)+1i*YY(i+1)-1i*YY(i);
    App.angle(i) = angle(ri);

    Boxr = [XX(i+1)-XX(i), YY(i+1)-YY(i)];
    b = [-YY(i+1)+YY(i), XX(i+1)-XX(i)]/abs(ri)*2*wk;
    rt = repmat([XX(i),YY(i)],5,1)+[b/2; b/2+Boxr; -b/2+Boxr; -b/2; b/2];
    App.box(:,:,i) = rt;

    line(rt(:,1)*zm,rt(:,2)*zm,'color','w','linew',1)

    
    if(i==1)
        roi(1) = min(rt(:,1));
        roi(2) = max(rt(:,1));
        roi(3) = min(rt(:,2));
        roi(4) = max(rt(:,2));
    else
        roi(1) = min([min(rt(:,1)),roi(1)]);
        roi(2) = max([max(rt(:,1)),roi(2)]);
        roi(3) = min([min(rt(:,2)),roi(3)]);
        roi(4) = max([max(rt(:,2)),roi(4)]);
    end

end

print(gcf,fullfile(resultpath,'Selected_ROI_2'),'-dpng','-r300');

h2 = figure('Name','Analysis ROI Enlarge','Position',[100,100,800,800]);
copyobj(ax,h2);
xlim([roi(1)-1 roi(2)+1]*zm)
ylim([roi(3)-1 roi(4)+1]*zm)
set(gcf,'color','w','Position',[850   100   (roi(2)-roi(1)+2)*30   (roi(4)-roi(3)+2)*30])
print(gcf,fullfile(resultpath,'Selected_ROI_Enlarge'),'-dpng','-r300');

%%  Z range selection

if ~isfield(para,'Lz')
    nref = round(length(XX)/2);
    xyzT = [cos(App.angle(nref)) sin(App.angle(nref)) 0; -sin(App.angle(nref)) cos(App.angle(nref)) 0; 0 0 1]*[recony'-hh;reconx'-hh;reconz'];

    reconxT = xyzT(2,:)'+hh;
    reconyT = xyzT(1,:)'+hh;
    reconzT = xyzT(3,:)';

    xyzTL = [cos(App.angle(nref)) sin(App.angle(nref)) 0; -sin(App.angle(nref)) cos(App.angle(nref)) 0; 0 0 1]*[XX'-hh;YY'-hh;zeros(size(XX'))];

    TxL = xyzTL(1,:)'+hh;
    TyL = xyzTL(2,:)'+hh;

    TxL_zm = TxL*zm;
    TyL_zm = TyL*zm;

    wkt = wk;

    Lxt = round(TyL_zm(nref)-wkt*zm : TyL_zm(nref)+wkt*zm);
    Lyt = TxL_zm(2)/zm : TxL_zm(length(XX)-2)/zm;
   masksub = reconyT>Lyt(1)&reconyT<Lyt(end) & reconxT>Lxt(1)/zm & reconxT<Lxt(end)/zm;

    Lzt = min(reconzT(masksub))/psz*zm-zm : max(reconzT(masksub))/psz*zm+zm;

    reconxTs = reconxT(masksub);
    reconzTs = reconzT(masksub);

    sz2 = ceil(max([(Lzt(end)-Lzt(1)) , (Lxt(end)-Lxt(1))])/zm);
    [rch1, gch1, bch1] = srhist_color_ca(sz2,zm,reconzTs/psz-mean(reconzTs/psz)+sz2/2,reconxTs-mean(reconxTs)+sz2/2,reconzTs,cm,ca0);

    rchsm_1 = gaussf(rch1,sigma );
    gchsm_1 = gaussf(gch1,sigma );
    bchsm_1 = gaussf(bch1,sigma);

    rchsmst_1 = imstretch_linear(rchsm_1,0,Recon_color_highb*10,0,255);
    gchsmst_1 = imstretch_linear(gchsm_1,0,Recon_color_highb*10,0,255);
    bchsmst_1 = imstretch_linear(bchsm_1,0,Recon_color_highb*10,0,255);


    colorim_1 = double(joinchannels('RGB',rchsmst_1,gchsmst_1,bchsmst_1));
        tmp = colorim_1(colorim_1>2);
        [val,ed] = histcounts(tmp,100,'Normalization','cdf');
        indct = find(val>0.99);
        normf = ed(indct(1));

    %%
        h3 = figure('Name','Select Z range','color','k','Position',[1200   180   700  700]);
        h3.InvertHardcopy = 'off';
        axes('Position',[0,0,1,1]);
    image(colorim_1./normf);
        axis equal 
        axis off

    [~,y2_zm] = myginput(2);
    line([Lxt(1)-mean(Lxt)+sz2*zm/2;Lxt(1)-mean(Lxt)+sz2*zm/2;Lxt(end)-mean(Lxt)+sz2*zm/2;Lxt(end)-mean(Lxt)+sz2*zm/2;Lxt(1)-mean(Lxt)+sz2*zm/2],[y2_zm;y2_zm(2);y2_zm(1);y2_zm(1)],'color','y','linew',2)
    drawnow
    Lz = round(mean(reconzTs/psz)*zm-zm*sz2/2+min(y2_zm)-2 : mean(reconzTs/psz)*zm-zm*sz2/2+max(y2_zm)+2);
else
    Lz = para.Lz;
end

%% Split and align the cross sections and save 

imS = zeros(size(rch,1),size(rch,2),3);

ROI = zeros(4,2);
Spp = struct;
Spp.reconxTs = [];
Spp.reconyTs = [];
Spp.reconzTs = [];
Spp.secN = [];
Spp.box = zeros(size(App.box));


pause(5)
close all

for ii =  1:length(XX)-1
    disp(['Crop section ' num2str(ii) ' / ' num2str(length(XX)-1)]);
    
    xyzTL = [cos(App.angle(ii)) sin(App.angle(ii)) 0; -sin(App.angle(ii)) cos(App.angle(ii)) 0; 0 0 1]*[App.box(:,:,ii)'-hh;zeros(1,5)];
    xyzTL(1:2,:) = xyzTL(1:2,:)+hh;
    cc = [(xyzTL(1,1)+xyzTL(1,4))/2, (xyzTL(2,1)+xyzTL(2,4))/2];
    

    xyzT = [cos(App.angle(ii)) sin(App.angle(ii)) 0; -sin(App.angle(ii)) cos(App.angle(ii)) 0; 0 0 1]*[recony'-hh;reconx'-hh;reconz'];
    xyzT(1:2,:) = xyzT(1:2,:)+hh;
    
    if(ii==1)
        rx = xyzTL(1,:)';
        ry = xyzTL(2,:)';
        reconxT = xyzT(2,:)';
        reconyT = xyzT(1,:)';
        reconzT = xyzT(3,:)';
        
        ROI(1,:) = [rx(1),ry(1)]*zm;
        ROI(4,:) = [rx(4),ry(4)]*zm;
        
    else
        rx = xyzTL(1,:)'-cc(1)+c0(1);
        ry = xyzTL(2,:)'-cc(2)+c0(2);
        reconxT = xyzT(2,:)'-cc(2)+c0(2);
        reconyT = xyzT(1,:)'-cc(1)+c0(1);
        reconzT = xyzT(3,:)';

    end
    
    ROI(2,:) = [rx(2),ry(2)]*zm;
    ROI(3,:) = [rx(3),ry(3)]*zm;

    c0 = [(rx(2)+rx(3))/2, (ry(2)+ry(3))/2];

    %%
    tt = 4;
    Ly2 = floor(rx(1)*zm) : ceil(rx(2)*zm);
    Lx2 = floor(ry(3)*zm) : ceil(ry(1)*zm);

    msk = reconyT>=(Ly2(1)-tt)/zm & reconyT<=(Ly2(end)+tt)/zm     & reconxT>=(Lx2(1)-tt)/zm & reconxT<=(Lx2(end)+tt)/zm     & reconzT>=Lz(1)/zm*psz & reconzT<=Lz(end)/zm*psz;
    
    reconxTs = reconxT(msk);
    reconyTs = reconyT(msk);
    reconzTs = reconzT(msk);
    
    n = length(Spp.reconxTs);
    Spp.reconxTs(n+1:n+length(reconxTs)) = reconxTs;
    Spp.reconyTs(n+1:n+length(reconxTs)) = reconyTs;
    Spp.reconzTs(n+1:n+length(reconxTs)) = reconzTs;
    Spp.secN(n+1:n+length(reconxTs)) = ii;
    Spp.box(:,:,ii) = [rx,ry];
    
    reconxTs_r = reconxTs - min(Lx2-tt)/zm;
    reconyTs_r = reconyTs - min(Ly2-tt)/zm;
    
    sz2 = ceil(max([max(Ly2)-min(Ly2)+2*tt,max(Lx2)-min(Lx2)+2*tt]/zm));
    
    [rch1, gch1, bch1] = srhist_color_ca(sz2,zm,reconxTs_r,reconyTs_r,reconzTs,cm,ca0);

    rchsm_1 = gaussf(rch1,sigma(1));
    gchsm_1 = gaussf(gch1,sigma(1));
    bchsm_1 = gaussf(bch1,sigma(1));

    rchsmst_1 = imstretch_linear(rchsm_1,0,Recon_color_highb/0.6,0,255);
    gchsmst_1 = imstretch_linear(gchsm_1,0,Recon_color_highb/0.6,0,255);
    bchsmst_1 = imstretch_linear(bchsm_1,0,Recon_color_highb/0.6,0,255);

    im_s = double(joinchannels('RGB',rchsmst_1,gchsmst_1,bchsmst_1));
    imS(Lx2,Ly2,:) = im_s(Lx2-min(Lx2-tt)+1,Ly2-min(Ly2-tt)+1,:);

end

%%
save(fullfile(resultpath,'Result_mid.mat'),'xy_zm','XX','YY','ca0','Lz','ROI','imS','App','Spp','para','subname')
clear mex
