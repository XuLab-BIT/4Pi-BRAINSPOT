function obj = genZKresult(obj)
% genZKresult - expand phase retrieved pupil function into zernike polynomials. 
%   It uses the object, 'obj.Z', from Zernike_Polynomials class.
%
%   see also Zernike_Polynomials
%% pupil a
pupil_phase=obj.PRstruct1.Pupil.phase;
pupil_mag=obj.PRstruct1.Pupil.mag;
n=obj.PRstruct.RefractiveIndex;
R=obj.PSFsize;
Z_N=obj.ZernikeorderN;
if isreal(pupil_phase)
    R_aber=pupil_phase.*obj.NA_constrain;
else
    R_aber=angle(pupil_phase).*obj.NA_constrain;
end
U=pupil_mag.*cos(R_aber).*obj.NA_constrain;
V=pupil_mag.*sin(R_aber).*obj.NA_constrain;
complex_Mag=U+1i*V;

[CN_complex,pupil_complexfit]=obj.Z.fitzernike(complex_Mag,'mag',   Z_N, R);
[CN_phase,pupil_phasefit]    =obj.Z.fitzernike(R_aber,     'phase', Z_N, R);
[CN_mag,pupil_magfit]        =obj.Z.fitzernike(pupil_mag,  'mag',   Z_N, R);

obj.PRstruct1.Zernike_phase=CN_phase;
obj.PRstruct1.Zernike_phaseinlambda=CN_phase./2./pi;
obj.PRstruct1.Zernike_mag=CN_mag;
obj.PRstruct1.Zernike_complex=CN_complex;
obj.PRstruct1.Fittedpupil.complex=pupil_complexfit;
obj.PRstruct1.Fittedpupil.phase=pupil_phasefit;
obj.PRstruct1.Fittedpupil.mag=pupil_magfit;

%% pupil b
pupil_phase=obj.PRstruct2.Pupil.phase;
pupil_mag=obj.PRstruct2.Pupil.mag;
n=obj.PRstruct.RefractiveIndex;
R=obj.PSFsize;
Z_N=obj.ZernikeorderN;
if isreal(pupil_phase)
    R_aber=pupil_phase.*obj.NA_constrain;
else
    R_aber=angle(pupil_phase).*obj.NA_constrain;
end
U=pupil_mag.*cos(R_aber).*obj.NA_constrain;
V=pupil_mag.*sin(R_aber).*obj.NA_constrain;
complex_Mag=U+1i*V;

[CN_complex,pupil_complexfit]=obj.Z.fitzernike(complex_Mag,'mag',   Z_N, R);
[CN_phase,pupil_phasefit]    =obj.Z.fitzernike(R_aber,     'phase', Z_N, R);
[CN_mag,pupil_magfit]        =obj.Z.fitzernike(pupil_mag,  'mag',   Z_N, R);


obj.PRstruct2.Zernike_phase=CN_phase;
obj.PRstruct2.Zernike_phaseinlambda=CN_phase./2./pi;
obj.PRstruct2.Zernike_mag=CN_mag;
obj.PRstruct2.Zernike_complex=CN_complex;
obj.PRstruct2.Fittedpupil.complex=pupil_complexfit;
obj.PRstruct2.Fittedpupil.phase=pupil_phasefit;
obj.PRstruct2.Fittedpupil.mag=pupil_magfit;