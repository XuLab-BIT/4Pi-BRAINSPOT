%% Refine the pupil function and estimate the aberration 


function [probj,empupil] = PR4Pi_fromAveZ_4channels_Zshift(ims_Zcal_ave_planes,index_record_Zplanes,empupil)


%% create object and set input properties of PRPSF class
disp('Generate pupil function and estimate aberration based on Phase retrieved method');
probj = [];
probj.PRstruct.NA = empupil.NA; 
probj.PRstruct.Lambda = empupil.Lambda;
probj.PRstruct.RefractiveIndex = empupil.nMed;
probj.Pixelsize = empupil.Pixelsize;  
probj.PSFsize = 128;
probj.SubroiSize = empupil.imsz;    
probj.OTFratioSize = 60;  
probj.ZernikeorderN = empupil.ZernikeorderN;  

probj.IterationNum = 60;  
probj.IterationNumK = 10; 

probj.Phasediff = empupil.Phasediff;                            % phase difference between s- and p-polarizations
probj.Iratio = empupil.Iratio;                                % transmission ratio between top and bottom emission path, from 0 to 1
probj.Phi0 = empupil.Phi0;                                 % cavity phase
probj.Zoffset = empupil.Zoffset;    
probj.ModulationDepth = empupil.ModulationDepth;

% probj.phase1 = empupil.phase1;
% probj.phase2 = empupil.phase2;

%% prepare the PSFs and its positions
display('Prepare PR data...');
Zpos_plane1 = empupil.Z_pos + empupil.zshift; 

display(['Z-position index: ' num2str(index_record_Zplanes')]);

Zpos_plane1_sel = Zpos_plane1(index_record_Zplanes==1);    
ims_Zcal_ave_planes_sel = ims_Zcal_ave_planes(:,:,index_record_Zplanes==1,:);


probj.Zindstart = 1; %index of position
probj.Zindend = size(Zpos_plane1_sel,2);
probj.Zindstep = 1;

%% generate PR result
probj.Zpos = Zpos_plane1_sel;
input_ZData = ims_Zcal_ave_planes_sel;
probj.BeadData =  input_ZData;    
probj.Beadcenter = [size(input_ZData,1)/2-1 size(input_ZData,2)/2-1];
probj.BeadXYshift = [0 0];

%% prepare initial parameters 
probj = precomputeParam(probj);
probj = datapreprocess(probj);

%% Carry out 4Pi phase retrieval 

probj = genMpsf(probj); 
probj = phaseretrieve4pi_4channels_extraDefocus_v3(probj);

probj = genZKresult(probj);


%% shift XYZ
if empupil.Zshift_mode
    for ii = 1 : 3
        
        probj.BeadData = input_ZData;
        
        C4_1 = probj.PRstruct1.Zernike_phase(4);
        est_1 = fminsearch(@(x)fitdefocus(probj,x,C4_1),[0.2,1]);% calculate defocus
        C4_2 = probj.PRstruct2.Zernike_phase(4);
        est_2 = fminsearch(@(x)fitdefocus(probj,x,C4_2),[0.2,1]);% calculate defocus
        
        zshift_tmp = (est_2(1)-est_1(1)) / 2;

        
        CXY_1 = probj.PRstruct1.Zernike_phase([2,3]);
        xyshift_1 = CXY_1*probj.PRstruct.Lambda/(2*pi*probj.Pixelsize*probj.PRstruct.NA);% calculate XY shift
        CXY_2 = probj.PRstruct2.Zernike_phase([2,3]);
        xyshift_2 = CXY_2*probj.PRstruct.Lambda/(2*pi*probj.Pixelsize*probj.PRstruct.NA);% calculate XY shift
        xyshift = (xyshift_1 + xyshift_2) / 2;
        
        empupil.zshift = empupil.zshift + zshift_tmp;
        probj.Zpos = probj.Zpos + zshift_tmp;
        probj.BeadXYshift = probj.BeadXYshift-xyshift;
        
        probj = datapreprocess(probj);
                
        % updated PR
        probj = genMpsf(probj);
        probj = phaseretrieve4pi_4channels_extraDefocus_v3(probj);

        probj = genZKresult(probj);
        

    end
    
end
    
%% 


probj = findOTFparam_fixedsigma(probj);  
% probj = findOTFparam_fixedsigma(probj,empupil.blur_sigma);   
% probj.PRstruct.SigmaX
if empupil.Zshift_mode == 2
    empupil.zshift = 0;
end

%% generate figures for phase retrieval results


display(['Z shift is: ' num2str(empupil.zshift)]);

display(['Updated Zernike(Z5-Z9) of pupil A: ' num2str(probj.PRstruct1.Zernike_phase(5:9))]);
display(['Updated Zernike(Z5-Z9) of pupil B: ' num2str(probj.PRstruct2.Zernike_phase(5:9))]);





