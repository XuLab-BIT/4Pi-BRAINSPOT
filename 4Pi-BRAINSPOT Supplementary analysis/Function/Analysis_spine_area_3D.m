function Analysis_spine_area_3D(para,subname,alpha)
%%
resultpath = [para.savepath, subname, '\'];

if (exist(fullfile(resultpath,'Result_mid.mat'),'file'))
    load(fullfile(resultpath,'Result_mid.mat'));
else
    error('Result_mid not found.');
end

fileID = fopen(fullfile(resultpath,'Area.txt'),'w');
fprintf(fileID,'%12s %18s %18s %18s %18s %18s\r\n','Position(nm)','long axis(nm)','short axis(nm)', 'mean axis (nm)','area (nm^2)', 'mean radius (nm)');

warning off
%%
%%  Parameter
s = 2;

k = para.k;


psz = para.psz;
zm = para.zm;
sigma = para.sigma;


Recon_color_highb = para.Recon_color_highb;

tk = para.tk/psz;      %% convert to pixel    


cm = para.cm;


savepath = para.savepath;
resultpath = [para.savepath, subname, '\'];
if ~exist(resultpath,'dir')
    mkdir(resultpath)
end

resultpath2 = [savepath, subname, '\Area\'];
if ~exist(resultpath2,'dir')
    mkdir(resultpath2)
end

TagS = 0;

sz2 = ceil(max([ROI(1,2)-ROI(3,2),max(Lz)-min(Lz)]/zm));
a1 = min([(max(Lz)-min(Lz))*s,sz2*zm*s]);
a2 = min([ceil((ROI(1,2)-ROI(3,2))*s)+1,sz2*zm*s]);

%%  Analysize section by section
Ar0 = NaN(1,length(XX)-1);
Dm = Ar0;
DM = Ar0;
D_hist = Ar0;
x0 = Ar0; 
z0 = Ar0;

for i3 = 1:length(XX)-1
%% Visualize XZ cross section

    disp(['Analyze section ' num2str(i3) ' / ' num2str(length(XX)-1)]);
    msksub = Spp.secN>=i3-k & Spp.secN<=i3+k;
    
    reconxTsW = Spp.reconxTs(msksub)';
    reconzTsW = Spp.reconzTs(msksub)';
    
    reconxTsW_r = reconxTsW - ROI(3,2)/zm;
    reconzTsW_r = (reconzTsW - Lz(1)/zm*psz)/psz ;
    reconyTsW = Spp.reconyTs(msksub)';
    
    
    [rch1W, gch1W, bch1W] = srhist_color_ca(sz2,zm*s,reconzTsW_r,reconxTsW_r,reconzTsW,cm,ca0);

    rchsm_1W = gaussf(rch1W,sigma*s);
    gchsm_1W = gaussf(gch1W,sigma*s);
    bchsm_1W = gaussf(bch1W,sigma*s);

    M = Recon_color_highb* max(rchsm_1W(:)+gchsm_1W(:)+bchsm_1W(:));
    rchsmst_1W = imstretch_linear(rchsm_1W,0,M,0,255);
    gchsmst_1W = imstretch_linear(gchsm_1W,0,M,0,255);
    bchsmst_1W = imstretch_linear(bchsm_1W,0,M,0,255);

    colorim_1W = double(joinchannels('RGB',rchsmst_1W,gchsmst_1W,bchsmst_1W));
        tmp = colorim_1W(colorim_1W>2);
        [val,ed] = histcounts(tmp,100,'Normalization','cdf');
        indct = find(val>0.995);
    if(~isempty(indct))
        normf = ed(indct(1));
        
        imlW = colorim_1W(1:a1, 1:a2, :)./normf;
        h3 = figure('visible','off','Name',['Slice_' num2str(i3)],'color','k','Position',[1000   180   400/(max(Lz)-min(Lz))*(ROI(1,2)-ROI(3,2))   400]);
        h3.InvertHardcopy = 'off';
        axes('Position',[0,0,1,1]);
        image(imlW);
        axis equal 
        axis off
        if (TagS==1)
            set(h3, 'visible', 'on'); 
        else
            imwrite(imlW,fullfile(resultpath2,['\Selected_YZ_C' num2str(i3) '.tiff']));
            close(h3)
        end
        
    end

    

    %% Remove outliers and initiate alphashape algorithm
    
    thc = 2;
    [Nx,edx] = histcounts(reconxTsW,min(reconxTsW):thc/zm:max(reconxTsW));
    [Nz,edz] = histcounts(reconzTsW,min(reconzTsW):thc*psz/zm:max(reconzTsW));
    
    th = min([(thc+1)*(k+1), floor(thc/100*sum(msksub)*(k+1))]);
    idx1 = find(Nx>=th,1,'first');
    idx2 = find(Nx>=th,1,'last');
    idz1 = find(Nz>=th,1,'first');
    idz2 = find(Nz>=th,1,'last');
    msk2 = reconxTsW>edx(idx1) & reconxTsW<edx(idx2) & reconzTsW>edz(idz1) & reconzTsW<edz(idz2);
    
    shp = alphaShape(double(reconxTsW(msk2)),double(reconzTsW(msk2)/psz),alpha);
   
    if (numRegions(shp)>0)
    %%  conduct alphashape algorithm for pattern recognition and determine pattern center
        msk4 = zeros(size(msk2));
        for i4 = 1:numRegions(shp)
            if i4==1
                dc1 = area(shp,i4);
                msk4 = inShape(shp,double(reconxTsW),double(reconzTsW/psz),i4);
                
            else
                dc = area(shp,i4);
                if (dc>=0.6*dc1 && dc<=1.4*dc1)
                    msk4 = msk4 | inShape(shp,double(reconxTsW),double(reconzTsW/psz),i4);
                elseif(dc>dc1)
                    dc1 = dc;
                    msk4 = inShape(shp,double(reconxTsW),double(reconzTsW/psz),i4);
                end
            end
        end
        
        
        shp2 = alphaShape(double(reconxTsW(msk4)),double(reconzTsW(msk4)/psz),2);
        
        h4 = figure('visible','off','Name',['Slice_' num2str(i3) '_ellipse'],'color','w','Position',[1000   650   400/(Lz(end)-Lz(1))*(ROI(1,2)-ROI(3,2))   400]);
        ax = gca;
        hold on
        set(ax, 'YDir','reverse')
        plot(shp2,'Edgecolor','none','FaceAlpha',0.6)
        daspect([1,1,1])
        box on
        xlim([ROI(3,2)/zm, ROI(1,2)/zm])
        ylim([Lz(1)/zm, Lz(end)/zm])
        ax.XTick = linspace(ROI(3,2)/zm+30/psz,ROI(1,2)/zm-30/psz,3);
        ax.YTick = linspace(Lz(1)/zm+30/psz,Lz(end)/zm-30/psz,3);
        xticklabels(round((ax.XTick-ax.XTick(2))*psz,0));
        yticklabels(round((ax.YTick-ax.YTick(2))*-psz,0));
        set(gca,'fontsize',25,'linewidth',3,'fontweight','bold')
        
        set( ax,'NextPlot','add' );
        plot(double(reconxTsW),double(reconzTsW/psz),'k.','Markersize',5)

        [~,P] = boundaryFacets(shp2);
        Ek = fit_ellipse([P(:,1);double(reconxTsW(msk4))],[P(:,2);double(reconzTsW(msk4)/psz)],ax);
        
        if TagS==1
            set(h4, 'visible', 'on'); 
        else
            print(h4,fullfile(resultpath2,['\Selected_YZ_S' num2str(i3)]),'-dpng','-r300');
            close(h4)
        end
        
        %%   Calculate geometric mean radius from center to localizations
        
        if ~isempty(Ek) && ~isempty(Ek.X0_in)
            x0(i3) = Ek.X0_in;
            z0(i3) = Ek.Y0_in;
            rr = sqrt((reconxTsW(msk2)-Ek.X0_in).^2 + (reconzTsW(msk2)/psz-Ek.Y0_in).^2)*psz;
            
            bin = sigma(1)*psz*sqrt(2*log(2));
            edr = double(0 : bin/zm : 320);     % nm
            Edr = edr(1:end-1)+diff(edr(1:2))/2;
            xf = linspace(min(Edr),max(Edr),20*numel(Edr));
    
            Nr = histcounts(rr,edr);
            options = fitoptions('gauss1','Lower',[0.8*max(Nr) mean(rr)-std(rr) 0],'Upper',[1.2*max(Nr) mean(rr)+std(rr) inf]);
            ff = fit(double(Edr'),double(Nr'),'gauss1',options);
            yf = feval(ff,xf)';
            b1 = ff.b1;
          
            h5 = figure('visible','off','Name','Line profile','color','w');
            hold on
            box on
            bar(Edr,Nr/max(Nr),1,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.4,'EdgeColor','none')
            plot(xf,yf/max(Nr),'color',[0.8500 0.3250 0.0980],'LineWidth',3)
            plot([b1 b1], [0 ff.a1/max(Nr)], 'k','LineWidth',3)

            yticks([0 1])
            xlabel('Position (nm)')
            ylabel('Norm. Count')
            xlim([min(Edr) max(Edr)])      
            ylim([0 1.2])
           xlim([min(xf) 80])
            set(gcf,'color','w','Position',[350   400   550 350])   
            set(gca,'position',[0.15 0.3 0.8 0.65],'fontsize',25,'linewidth',2.5,'fontweight','bold')
            box off
            
            if TagS==1
                set(h5, 'visible', 'on'); 
            else
                print(h5,fullfile(resultpath2,['\Radius_YZ_S' num2str(i3)]),'-dpng','-r300');
                close(h5)    
            end
            %%
            Dm(i3) = Ek.short_axis*psz/2;
            DM(i3) = Ek.long_axis*psz/2;
            D_hist(i3) = ff.b1;
            Ar0(i3) = ff.b1.^2*pi;
            if TagS==0
                fprintf(fileID,'%12.2f %18.2f %18.2f %18.2f %18.2f %18.2f\r\n', (i3-0.5)*para.tk, DM(i3), Dm(i3), sqrt(DM(i3)*Dm(i3)), Ar0(i3), D_hist(i3));
            end
        else
            if TagS==0
                fprintf(fileID,'%12.2f \r\n', (i3-0.5)*para.tk);
            end
        end


    end
    %}
    drawnow
    
end
fclose(fileID);

%%  remove outlier point in volume estimation
x = (1:length(Ar0))*tk*psz;
Ar = Ar0;

figure
plot(1:length(Ar0),Ar0)
ylim([0 inf])

outlier = input('Input outlier position array (ex: [3,29])/(Empty for no) :  ');
if ~isempty(outlier)
    Ar(outlier) = nan;
end

jj = ~isnan(Ar);

f = fit(x(jj)',Ar(jj)','smoothingspline','SmoothingParam',1e-5);
xf = linspace(min(x),max(x),numel(x)*10);
yf = feval(f,xf);
tl = find(yf>=0,1,'last');
[~,ml] = max(yf(round(tl/2):tl(end)));
ml = ml+round(tl/2)-1;
yM = max([15,max(Ar(jj)/10^4)]);

figure
hold on
plot(xf,yf/10^4,'color',[0.9290 0.6940 0.1250],'linew',4)
plot(x,Ar/10^4,'.','color',[0 0.4470 0.7410],'linew',1,'Markersize',15)
xlim([min(xf) max(xf)])
xlabel('Position (nm)')
ylabel('Area (10^4 nm^2)')
xticks(round(linspace(ceil(min(xf)/100), floor(max(xf)/100-0.5), 3))*100)
legend({'Area','Fitting'},'Location','SouthEast')
set(gcf,'color','w','Position',[50   300   780 400])
set(gca,'fontsize',22,'linewidth',2)
print(gcf,fullfile(resultpath,'Selected_ROI_straight_Area2'),'-dpng','-r300');

tmp = imS(imS>2);
    [val,ed] = histcounts(tmp,100,'Normalization','cdf');
    indct = find(val>0.99);
    normf = ed(indct(1));
    h4 = figure('Name','Straighten pattern','color','k','Position',[50   180   300/(ROI(2,2)-ROI(3,2))*(ROI(2,1)-ROI(1,1))   300]);
    h4.InvertHardcopy = 'off';
    axes('Position',[0,0,1,1]);
image(imS./normf);
    axis equal 
    axis off
xlim(round([ROI(1,1)-1,ROI(2,1)+1]))
ylim(round([ROI(3,2)-1,ROI(2,2)+1]))
print(h4,fullfile(resultpath,'Selected_ROI_straight_XYs'),'-dpng','-r300');

for ii = 1:length(XX)-1
    figure(h4)
    line(Spp.box(:,1,ii)*zm,Spp.box(:,2,ii)*zm,'color','w','linew',2)
end
print(h4,fullfile(resultpath,'Selected_ROI_straight_XYs2'),'-dpng','-r300');


%%  Save
save(fullfile(resultpath,'Result_final.mat'),'xy_zm','ca0','XX','YY','alpha','Lz','ROI','imS','x','x0','z0','Ar0','Ar','Spp','para','subname')

ssOL = num2str(find(isnan(Ar)),' %d, ');
ssOL = ssOL(1:end-1);

fileID = fopen(fullfile(resultpath,'Parameter.txt'),'w');
fprintf(fileID,'%s  %s\r\n', 'Filepath :', para.filepath);
fprintf(fileID,'%s  %2.2f \r\n','Pixel size (nm) :', psz);
fprintf(fileID,'%s  %2.2f \r\n','Enlarge ratio :', zm);
fprintf(fileID,'%s  %2.2f \r\n','Section thickness (nm) :', para.tk*(2*k+1));
fprintf(fileID,'%s  %2.2f \r\n','Section width (nm) :', 2*para.wk);
fprintf(fileID,'%s  %2.2f \r\n','Running window shift (nm) :', para.tk);
fprintf(fileID,'%s  %2.2f \r\n','Alpha value :', alpha);
fprintf(fileID,'%s  %2.2f \r\n','Cross section rendered width :', a2*para.psz/para.zm/s);
fprintf(fileID,'%s  %2.2f \r\n','Cross section rendered height :', a1*para.psz/para.zm/s);
fprintf(fileID,'%s  %2.0f \r\n','Z select Upper bound :', max(Lz));
fprintf(fileID,'%s  %2.0f \r\n','Z select Lower bound :', min(Lz));
fprintf(fileID,'\r\n');
fprintf(fileID,'Estimation based on smoothingspline fitting result : \r\n');
fprintf(fileID,'%s  %s \r\n','Outlier point :', ssOL);
fprintf(fileID,'%s  %2.2f \r\n','Total length (nm) :', xf(tl));
fprintf(fileID,'%s  %2.2f \r\n','Total volume (nm^3) :', trapz(xf(1:tl),yf(1:tl)));
fprintf(fileID,'%s  %2.2f \r\n','Maximal area of head (nm^2) :', yf(ml));
fprintf(fileID,'%s  %2.2f \r\n','Minimal area of neck (nm^2) :', min(yf(1:ml)));
fclose(fileID);

% pause(2)