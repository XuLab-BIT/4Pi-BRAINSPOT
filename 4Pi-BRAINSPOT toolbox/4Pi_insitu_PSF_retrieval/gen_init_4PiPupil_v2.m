%Generate initial pupil from zernike value

function probj = gen_init_4PiPupil_v2(empupil,label_top,label_bot)

disp('generate initial pupil');

% Initialize
R = 128;
PRstruct = [];
PRstruct.NA = empupil.NA;                   % numerical aperture of the objective lens
PRstruct.Lambda = empupil.Lambda;           % center wavelength of the emission band spass filter, unit is micron
PRstruct.RefractiveIndex = empupil.nMed;    % refractive index of immersion oil
PRstruct.Pupil.phase = zeros(R,R);         
PRstruct.Pupil.mag = zeros(R,R);
PRstruct.SigmaX = empupil.blur_sigma;       
PRstruct.SigmaY = empupil.blur_sigma;
PRstruct1 = PRstruct;                             
PRstruct2 = PRstruct;                               

% Input magnitude and phase
magZ = zeros(1,empupil.Zernike_sz);
magZ(1) = 1;
phaseZ_1 = zeros(1,empupil.Zernike_sz);
phaseZ_1([5:25]) = label_top;  
phaseZ_2 = zeros(1,empupil.Zernike_sz);
phaseZ_2([5:25]) = label_bot;  

PRstruct1.Zernike_phase = phaseZ_1;
PRstruct1.Zernike_mag = magZ;
PRstruct2.Zernike_phase = phaseZ_2;
PRstruct2.Zernike_mag = magZ;


% generate pupil
probj = PSF_4pi(PRstruct1);                       
probj.Pixelsize = empupil.Pixelsize;               % pixel size on the sample plane, unit is micron
probj.PSFsize = R;                                 % image size used for PSF generation
probj.nMed = empupil.nMed;

probj.Phasediff = empupil.Phasediff;                      % phase difference between s- and p-polarizations
probj.Iratio = empupil.Iratio;                                       % transmission ratio between top and bottom emission path
probj.ModulationDepth = empupil.ModulationDepth;                      % modulation strength of interferometric PSFs
probj.Phi0 = empupil.Phi0;                                      % cavity phase
probj.Zoffset = empupil.Zoffset;

probj.gen2Pupil(PRstruct1,PRstruct2);             


probj.PRstruct1.Pupil.phase = exp(1i.*probj.Pupila.phase);
probj.PRstruct1.Pupil.mag = probj.Pupila.mag;
probj.PRstruct2.Pupil.phase = exp(1i.*probj.Pupilb.phase);
probj.PRstruct2.Pupil.mag = probj.Pupilb.mag;


