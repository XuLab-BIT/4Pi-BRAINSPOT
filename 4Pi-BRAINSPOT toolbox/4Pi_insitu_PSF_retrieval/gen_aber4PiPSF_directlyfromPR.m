
%% Generate reference Z-postions' PSFs in 4 channels from pupil function
%  Edited by Fan Xu


function [ref_planes,empupil] = gen_aber4PiPSF_directlyfromPR(probj,empupil,label_top,label_bot)


%% reference Z image in plane1 and plane2

disp('Generate reference PSFs in different Z positions of two channels from pupil function');
imsz = empupil.imsz;  
Numimage = size(empupil.Z_pos,2);   
ref_planes = zeros(imsz,imsz,Numimage,4);


num = 1;
for ii = 1 : Numimage
    
    %% calculate real Z image in plane1 and plane2
    Ztrue_plane1 = empupil.Z_pos(ii) + empupil.zshift;
        
    %% get pupil from Phase retrieval
    PRstruct = probj.PRstruct;
    pxsz = probj.Pixelsize;
    R = 128;
    PRstruct.SigmaX = empupil.blur_sigma;
    PRstruct.SigmaY = empupil.blur_sigma;
    PRstruct.Pupil.phase = zeros(R,R);
    PRstruct.Pupil.mag = zeros(R,R);
    
   
    %% generate reference Z images in Plane1
    psfobj = PSF_4pi(PRstruct);
    psfobj.Xpos = 0;
    psfobj.Ypos = 0;
    psfobj.Zpos = Ztrue_plane1;    

    psfobj.Boxsize = imsz;
    psfobj.Pixelsize = pxsz; % micron
    psfobj.PSFsize = R;
    psfobj.nMed = probj.PRstruct.RefractiveIndex;
    
    psfobj.Phasediff = empupil.Phasediff;                        % phase difference between s- and p-polarizations
    psfobj.Iratio = empupil.Iratio;                                % transmission ratio between top and bottom emission path
    psfobj.ModulationDepth = empupil.ModulationDepth;                       % modulation strength of interferometric PSFs
    psfobj.Phi0 = empupil.Phi0;                              % cavity phase
    psfobj.Zoffset = empupil.Zoffset;
    
    

    PRstruct1 = probj.PRstruct1;
    PRstruct2 = probj.PRstruct2;
    
    psfobj.precomputeParam();
    psfobj.set2Pupil(PRstruct1,PRstruct2);
    psfobj.setUnifMag();  
    psfobj.genPupil_4pi_2();
    psfobj.genPSF_4pi_md()

        
    %% add I, bg
    I = permute(repmat(5000, [1, imsz, imsz]), [2,3,1]);
    bg = 0;

    label_ch = {'p1','s2','p2','s1'};
    for nn = 1:4
        psfobj.PSFs = psfobj.PSF4pi.(label_ch{nn});
        psfobj.scalePSF();
        psf = psfobj.ScaledPSFs;
        img = psf.*I;
        
        ref_planes(:,:,num,nn) = sum(img, 3)+ bg;
                
    end
        
    num = num + 1;
    
end

empupil.phase1 = PRstruct1.Pupil.phase;
empupil.phase2 = PRstruct2.Pupil.phase;


for nn = 1 : 4
    for ii = 1 : size(ref_planes,3)
        tmp_ch = ref_planes(:,:,ii,nn);
        tmp_ch = reshape(zscore(tmp_ch(:)),size(tmp_ch,1),size(tmp_ch,2));
        ref_planes(:,:,ii,nn) = tmp_ch; 
    end
end



