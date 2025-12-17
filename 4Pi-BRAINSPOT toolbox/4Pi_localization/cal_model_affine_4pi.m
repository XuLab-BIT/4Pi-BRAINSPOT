% 4 channel-specific model in 4pi configuration
function [tform_fitting_all, tMatrix_noTranslation] = cal_model_affine_4pi(tform_all,ref_point)

num_tform = length(tform_all);
tform_fitting_all = [];
tMatrix_noTranslation = [];
for ii = 2 : num_tform

    invertform = invert(tform_all{ii});
    
    [x,y] = transformPointsForward(invertform,ref_point(1),ref_point(2));
    
    x_offset = ref_point(1)-x;
    y_offset = ref_point(2)-y;
    
    tform_fitting = invertform;
    tform_fitting.T(3,1) = tform_fitting.T(3,1) + x_offset;
    tform_fitting.T(3,2) = tform_fitting.T(3,2) + y_offset;
    
    tform_fitting_all{ii-1} = tform_fitting;
    
    tMatrix_tmp = single(zeros(3,2));
    tMatrix_tmp(1,1) = tform_fitting.T(1,1);
    tMatrix_tmp(1,2) = tform_fitting.T(1,2);
    tMatrix_tmp(2,1) = tform_fitting.T(2,1);
    tMatrix_tmp(2,2) = tform_fitting.T(2,2);
    tMatrix_tmp(3,1) = tform_fitting.T(3,1);
    tMatrix_tmp(3,2) = tform_fitting.T(3,2);
    
    tMatrix_noTranslation = cat(2,tMatrix_noTranslation,tMatrix_tmp);
end
