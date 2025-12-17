%% Generate pupil function based on 4Pi phase retrieval

%%
function obj = phaseretrieve4pi_4channels_extraDefocus_v3(obj)

phib = obj.Phasediff;   % phase difference between s- and p-polarizations
Iratio = obj.Iratio;       % transmission ratio between top and bottom emission path
phi0 = obj.Phi0;        % cavity phase
Mdepth = obj.ModulationDepth;

Mpsf_select = obj.Mpsf_extend;
N_ch = obj.DatadimCh;

z = obj.Zpos;
zind = [obj.Zindstart:obj.Zindstep:obj.Zindend];
N = length(zind);
n = obj.PRstruct.RefractiveIndex;
Freq_max = obj.PRstruct.NA/obj.PRstruct.Lambda;
NA_constrain = obj.k_r<Freq_max;
k_z = sqrt((n/obj.PRstruct.Lambda)^2-obj.k_r.^2).*NA_constrain;
Fig = NA_constrain;

pupil_mag = Fig/sum(sum(Fig));

R = obj.PSFsize;


%% Initial phase
if isfield(obj,'phase1') && isfield(obj,'phase2')
    pupilA_phase = obj.phase1;
    pupilB_phase = obj.phase2;
else    
    pupilA_phase = ones(R,R);
    pupilB_phase = ones(R,R);
end

pupilA_mag = pupil_mag;
pupilB_mag = pupil_mag;
%%

for k = 1:obj.IterationNum

    N_pairs = N;
    RpupilA_mag = zeros(R,R,N_pairs);
    RpupilB_mag = zeros(R,R,N_pairs);
    RpupilA_phase = zeros(R,R,N_pairs);
    RpupilB_phase = zeros(R,R,N_pairs);
    for i = 1:1:N_pairs 
        Rpupil = zeros(R,R,N_ch);
        for j = 1:N_ch 
            
            defocusphaseA = exp(-2.*pi.*1i.*(z(zind(i))+obj.Zoffset).*k_z);%top
            defocusphaseB = exp(2.*pi.*1i.*z(zind(i)).*k_z); %bottom
            pupila = pupilA_mag.*pupilA_phase; %top
            pupilb = pupilB_mag.*pupilB_phase.*exp(1i.*phi0).*Iratio; %bottom 

            % updated 4pi pupil at differnt z positions
            if j == 1 %p1 channel
                pupil_complex = pupila.*exp(pi.*1i).*defocusphaseA + pupilb.*defocusphaseB;         
            elseif j == 2   %s2 channel
                pupil_complex = pupila.*defocusphaseA + pupilb.*defocusphaseB.*exp(1i*phib);
            elseif j == 3 %p2 channel
                pupil_complex = pupila.*defocusphaseA + pupilb.*defocusphaseB;
            else %s1 channel
                pupil_complex = pupila.*defocusphaseA.*exp(1i*pi) + pupilb.*defocusphaseB.*exp(1i*phib);
            end
            
            % incoherence part
            psfa = abs(fftshift(fft2(pupila.*defocusphaseA))).^2;  
            psfb = abs(fftshift(fft2(pupilb.*defocusphaseB))).^2; 

            % combine coherence and incoherence parts
            Fig1 = abs(fftshift(fft2(pupil_complex))).^2;       
            Fig1_all = Mdepth.*Fig1 + (1-Mdepth).*(psfa+psfb); 

            coeff_indx = sum(sum(abs(pupil_complex)));
            
            PSF0 = Fig1_all./sum(sum(Fig1_all));  
            Mpsfo = squeeze(Mpsf_select(:,:,zind(i),j));  
        
            
            Mpsfo = (Mpsfo - (1-Mdepth).*(psfa+psfb)./sum(sum(Fig1_all)))./Mdepth; %only consider coherence part
            Mpsfo(Mpsfo<=0) = 0;         
            Mpsfo = Mpsfo./sum(sum(Mpsfo));
            
            PSF1 = Fig1./sum(sum(Fig1));      
            if k>obj.IterationNumK
                Mask = (Mpsfo==0);
                Mpsfo(Mask) = PSF1(Mask);
            end    

            Rpsf = fft2(pupil_complex);  %FFT each complex section
            Rpsf_phase = Rpsf./abs(Rpsf);
            Fig2 = fftshift(sqrt(abs(Mpsfo)));  
            Mpsf = Fig2./sum(sum(Fig2));

            Rpupil(:,:,j) = ifft2((Mpsf.*Rpsf_phase));
            Rpupil(:,:,j) = Rpupil(:,:,j)./sum(sum(abs(Rpupil(:,:,j).*NA_constrain))).*coeff_indx; 

        end
        %% sovle(pupilA and pupilB)

        Rpupil_input = reshape(Rpupil,R*R,N_ch);
        
        AA = [exp(1i.*pi) Iratio.*exp(1i.*phi0);...
            1 Iratio.*exp(1i.*phib).*exp(1i.*phi0);...
            1 Iratio.*exp(1i.*phi0);...
            exp(pi.*1i) Iratio.*exp(1i.*phib).*exp(1i.*phi0)];
        

        Rpupil_ouput = AA \ Rpupil_input.';        
        Rpupil_ouput = reshape(Rpupil_ouput.',R,R,2); 
        Rpupil_ouput(Rpupil_ouput == 0) = 1e-8; 

        tmpA = Rpupil_ouput(:,:,1);
        tmpB = Rpupil_ouput(:,:,2);

        tmpA = tmpA./ sum(sum((abs(tmpA.*NA_constrain))));
        tmpB = tmpB./ sum(sum((abs(tmpB.*NA_constrain))));
        
        
        RpupilA_mag(:,:,i)=abs(tmpA);
        RpupilA_phase(:,:,i)=tmpA./RpupilA_mag(:,:,i).*conj(defocusphaseA);  %refocus
        
        RpupilB_mag(:,:,i)=abs(tmpB);
        RpupilB_phase(:,:,i)=tmpB./RpupilB_mag(:,:,i).*conj(defocusphaseB);


      
        
    end    
    %% average the pupils from pairs of PSFs 
        

    Fig3A = mean(RpupilA_phase.*RpupilA_mag,3);
    Fig3B = mean(RpupilB_phase.*RpupilB_mag,3);
    

    pupilA_phase = Fig3A./abs(Fig3A);
    pupilB_phase = Fig3B./abs(Fig3B);

    Fig5A = abs(Fig3A).*NA_constrain;
    Fig5B = abs(Fig3B).*NA_constrain;

    
    pupilA_mag = Fig5A;
    pupilA_mag = pupilA_mag./sum(pupilA_mag(:));

    pupilB_mag = Fig5B;
    pupilB_mag = pupilB_mag./sum(pupilB_mag(:));
    
    
    

end

% normalization
pupilA_mag = pupilA_mag.^2;
pupilA_mag = pupilA_mag./sum(sum(pupilA_mag));
pupilA_mag = sqrt(pupilA_mag); 

pupilB_mag = pupilB_mag.^2;
pupilB_mag = pupilB_mag./sum(sum(pupilB_mag));
pupilB_mag = sqrt(pupilB_mag); 

% generate phase retrieved PSF
psf_s1 = zeros(R,R,numel(z));
psf_s2 = zeros(R,R,numel(z));
psf_p1 = zeros(R,R,numel(z));
psf_p2 = zeros(R,R,numel(z));
psfa = zeros(R,R,numel(z));
psfb = zeros(R,R,numel(z));



for j=1:numel(z)
    defocusphaseA = exp(-2.*pi.*1i.*(z(zind(j))+obj.Zoffset).*k_z);%top
    defocusphaseB = exp(2.*pi.*1i.*z(zind(j)).*k_z); %bottom
    pupila = pupilA_mag.*pupilA_phase; %top
    pupilb = pupilB_mag.*pupilB_phase.*exp(1i.*phi0).*Iratio; %bottom
    
    pupila_complex = pupila.*defocusphaseA;
    tmpa = abs(fftshift(fft2(pupila_complex))).^2;
    psfa(:,:,j) = tmpa./R^2; % normalized PSF
    
    pupilb_complex = pupilb.*defocusphaseB;
    tmpb = abs(fftshift(fft2(pupilb_complex))).^2;
    psfb(:,:,j)=tmpb./R^2; % normalized PSF
   
    pupil_complex_s1 = pupila.*defocusphaseA.*exp(1i*pi) + pupilb.*defocusphaseB.*exp(1i*phib);
    tmp = abs(fftshift(fft2(pupil_complex_s1))).^2;
    psf_s1(:,:,j)=Mdepth.*tmp./R^2./4 + (1-Mdepth).*(psfa(:,:,j)+psfb(:,:,j))./4; % s1

    pupil_complex_s2 = pupila.*defocusphaseA + pupilb.*defocusphaseB.*exp(1i*phib);
    tmp2 = abs(fftshift(fft2(pupil_complex_s2))).^2;
    psf_s2(:,:,j)=Mdepth.*tmp2./R^2./4 + (1-Mdepth).*(psfa(:,:,j)+psfb(:,:,j))./4; % s2
    
    pupil_complex_p1 = pupila.*exp(pi.*1i).*defocusphaseA + pupilb.*defocusphaseB;
    tmp = abs(fftshift(fft2(pupil_complex_p1))).^2;
    psf_p1(:,:,j)=Mdepth.*tmp./R^2./4 + (1-Mdepth).*(psfa(:,:,j)+psfb(:,:,j))./4; % p1
    
    pupil_complex_p2 = pupila.*defocusphaseA + pupilb.*defocusphaseB;
    tmp2 = abs(fftshift(fft2(pupil_complex_p2))).^2;
    psf_p2(:,:,j)=Mdepth.*tmp2./R^2./4 + (1-Mdepth).*(psfa(:,:,j)+psfb(:,:,j))./4; % p2

end


% save pupil function and PSF in PRstruct
obj.PRstruct1.Pupil.phase = pupilA_phase;
obj.PRstruct1.Pupil.mag = pupilA_mag;
obj.PRstruct2.Pupil.phase = pupilB_phase;
obj.PRstruct2.Pupil.mag = pupilB_mag;
obj.PSFstruct.PRpsf = cat(4,psf_p1,psf_s2,psf_p2,psf_s1);

