% (C) Copyright 2022                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu

function genPupilfigs_4ch(obj,ImgType,datapath)
% genPRfigs - generate figures of phase retrieval result, including various PSFs, pupil 
% functions, and plots of zernike coefficients. 
%
%   Input parameter: ImgType - type of image to be generated, select from
%   'PSF', 'pupil' and 'zernike'
switch ImgType
    case 'PSF'        
        z = obj.Zpos;
        zind = [obj.Zindstart:obj.Zindstep:obj.Zindend];
        zind = zind(1:1:end); 
        N_ch = obj.DatadimCh;
        Ri = 24;    
        Ro = size(obj.Mpsf_extend,1);
        L = length(zind);
        label = {'p1','s2','p2','s1'};
        
        % measured PSFs
        min_color = min(obj.Mpsf_extend(:)); max_color = max(obj.Mpsf_extend(:));
        h1 = figure('Name','measured PSF at sampled z positions','position',[100,100,100*L,102*N_ch]);
        for ii = 1:N_ch
            Mpsf = squeeze(obj.Mpsf_extend(:,:,:,ii));  %measured PSF

            for jj = 1:L
                ha = axes('position',[(jj-1)/L,(N_ch-ii)/N_ch,1/L,1/N_ch],'parent',h1);
                imagesc(Mpsf(Ro/2-Ri/2+1:Ro/2+Ri/2,Ro/2-Ri/2+1:Ro/2+Ri/2,zind(jj)));
%                 imagesc(Mpsf(Ro/2-Ri/2+1:Ro/2+Ri/2,Ro/2-Ri/2+1:Ro/2+Ri/2,zind(jj)),[min_color max_color]);
                axis equal;axis off;
                if ii == 1
                    text(2,3, ['z=',num2str(z(zind(jj)),3),'\mum'],'color',[1,1,1],'fontsize',12);
                end
                if jj == 1
                    text(2,Ri-3, label{ii},'color',[1,1,1],'fontsize',12);
                end
            end
        end
        colormap(gray)
        saveas(gca,fullfile(datapath,'Assigned PSFs.tif'));   %edited by Fan Xu

        % PR PSFs
        min_color = min(obj.PSFstruct.PRpsf(:)); max_color = max(obj.PSFstruct.PRpsf(:));
        h1 = figure('Name','4Pi phase retrieved PSF at sampled z positions','position',[100,100,100*L,102*N_ch]);
        for ii = 1:N_ch
            Modpsf = squeeze(obj.PSFstruct.Modpsf(:,:,:,ii));   % Modpsf PR PSF  %obj.PSFstruct.PRpsf
            for jj = 1:L
                ha = axes('position',[(jj-1)/L,(N_ch-ii)/N_ch,1/L,1/N_ch],'parent',h1);
                imagesc(Modpsf(Ro/2-Ri/2+1:Ro/2+Ri/2,Ro/2-Ri/2+1:Ro/2+Ri/2,zind(jj)));
%                 imagesc(Modpsf(Ro/2-Ri/2+1:Ro/2+Ri/2,Ro/2-Ri/2+1:Ro/2+Ri/2,zind(jj)),[min_color max_color]);
                axis equal;axis off;
                if ii == 1
                    text(2,3, ['z=',num2str(z(zind(jj)),3),'\mum'],'color',[1,1,1],'fontsize',12);
                end
                if jj == 1
                    text(2,Ri-3, label{ii},'color',[1,1,1],'fontsize',12);
                end
            end
        end
        colormap(gray)
        saveas(gca,fullfile(datapath,'phase retrieved PSFs.tif'));   %edited by Fan Xu
        
    case 'pupil'
        % Top pupil function
        figure('Color',[1,1,1],'Name',' phase retrieved and Zernike fitted pupil function','Resize','on','Units','normalized','Position',[0.3,0.3,0.22,0.42])
        h1=[];
        RC=64;
        Rsub=63;
        h1(1)=subplot('Position',[0,0.5,1/2,1/2]);
        image(double(obj.PRstruct1.Pupil.mag(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h1(1))
        text(3,8,['PR pupil mag A'],'color',[1,1,1]);
        h1(2)=subplot('Position',[0,0,1/2,1/2]);
        image(double(obj.PRstruct1.Fittedpupil.mag(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h1(2))
        text(3,8,['Zernike pupil mag A'],'color',[1,1,1]);
        h1(3)=subplot('Position',[0.5,0.5,1/2,1/2]);
        tmp=angle(obj.PRstruct1.Pupil.phase);
        mag=obj.PRstruct1.Pupil.mag;
        mag(mag>0)=1;
        PRphase=tmp.*mag;

        image(double(PRphase(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h1(3))
        text(3,8,['PR pupil phase A'],'color',[0,0,0]);
        h1(4)=subplot('Position',[0.5,0,1/2,1/2]);
        image(double(obj.PRstruct1.Fittedpupil.phase(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h1(4))
        text(3,8,['Zernike pupil phase A'],'color',[0,0,0]);
        colormap(gray)
        axis(h1,'equal')
        axis(h1,'off')

        saveas(gca,fullfile(datapath,'phase retrieved and Zernike fitted pupil function (Top).tif'));   %edited by Fan Xu

        % Bottom pupil function 
        figure('Color',[1,1,1],'Name',' phase retrieved and Zernike fitted pupil function','Resize','on','Units','normalized','Position',[0.3,0.3,0.22,0.42])
        h2=[];
        RC=64;
        Rsub=63;
        h2(1)=subplot('Position',[0,0.5,1/2,1/2]);
        image(double(obj.PRstruct2.Pupil.mag(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h2(1))
        text(3,8,['PR pupil mag B'],'color',[1,1,1]);
        h2(2)=subplot('Position',[0,0,1/2,1/2]);
        image(double(obj.PRstruct2.Fittedpupil.mag(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h2(2))
        text(3,8,['Zernike pupil mag B'],'color',[1,1,1]);
        h2(3)=subplot('Position',[0.5,0.5,1/2,1/2]);
        tmp=angle(obj.PRstruct2.Pupil.phase);
        mag=obj.PRstruct2.Pupil.mag;
        mag(mag>0)=1;
        PRphase=tmp.*mag;

        image(double(PRphase(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h2(3))
        text(3,8,['PR pupil phase B'],'color',[0,0,0]);
        h2(4)=subplot('Position',[0.5,0,1/2,1/2]);
        image(double(obj.PRstruct2.Fittedpupil.phase(RC-Rsub:RC+Rsub,RC-Rsub:RC+Rsub)),'CDataMapping','scaled','Parent',h2(4))
        text(3,8,['Zernike pupil phase B'],'color',[0,0,0]);
        colormap(gray)
        axis(h2,'equal')
        axis(h2,'off')
        
        saveas(gca,fullfile(datapath,'phase retrieved and Zernike fitted pupil function (bottom).tif'));   %edited by Fan Xu
    case 'zernike'
        PlotZernikeC(obj.PRstruct1.Zernike_phase,'phase A');   % pupil a
        saveas(gca,fullfile(datapath,'Top phase_Zernike coefficinet.tif'));   %edited by Fan Xu

        PlotZernikeC(obj.PRstruct2.Zernike_phase,'phase B');   % pupil b
        saveas(gca,fullfile(datapath,'Bottom phase_Zernike coefficinet.tif'));   %edited by Fan Xu


end
end

function PlotZernikeC(CN_phase,type)
nZ=length(CN_phase);
vec=linspace(max(CN_phase)-0.1,min(CN_phase)+0.1,8);
dinv=vec(1)-vec(2);
ftsz=12;
figure('position',[200,200,700,300])
plot(CN_phase,'o-')
text(nZ+5,vec(1),['x shift: ', num2str(CN_phase(2),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(2),['y shift: ', num2str(CN_phase(3),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(3),['z shift: ', num2str(CN_phase(4),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(4),['Astigmatism: ', num2str(CN_phase(5),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(5),['Astigmatism(45^o): ', num2str(CN_phase(6),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(6),['Coma(x): ', num2str(CN_phase(7),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(7),['Coma(y): ', num2str(CN_phase(8),'%.2f')],'fontsize',ftsz);
text(nZ+5,vec(8),['Spherical: ', num2str(CN_phase(9),'%.2f')],'fontsize',ftsz);
xlim([0,nZ+30])
ylim([min(CN_phase)-dinv/2,max(CN_phase)+dinv/2])
set(gca,'fontsize',ftsz)
xlabel('Zernike coefficient number','fontsize',ftsz)
ylabel('Value','fontsize',ftsz)
title(type)


end
