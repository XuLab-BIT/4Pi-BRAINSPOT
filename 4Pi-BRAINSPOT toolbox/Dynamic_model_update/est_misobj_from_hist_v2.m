%% estimate the objective misalign parameter from measured histgraph

function [misobj_guess_keep,misobj_fit,totmask] = est_misobj_from_hist_v2(misobj_guess,z_guess,mcc_val)

%% rejection
ccmask = mcc_val < 0.6; 
zmask = z_guess > 0.5 | z_guess < -0.5; 
mismask = misobj_guess > 1.7 | misobj_guess < -1.7;
totmask = ccmask | zmask | mismask;
misobj_guess_keep = misobj_guess(~totmask);

%% parameter estimation

mean_misobj = mean(misobj_guess_keep);
misobj_input = misobj_guess_keep;

range_cal = 1;

tmp_misobj_low =  mean_misobj - range_cal;
tmp_misobj_high = mean_misobj + range_cal;

index_sel = misobj_input > tmp_misobj_low & misobj_input < tmp_misobj_high;
misobj_input_sel = misobj_input(index_sel);

pd = fitdist(misobj_input_sel,'normal');
tmp_mu = pd.mu;

for ii = 1 : 5
    tmp_misobj_low =  tmp_mu - range_cal;
    tmp_misobj_high = tmp_mu + range_cal;
    
    
    index_sel_2 = misobj_input > tmp_misobj_low & misobj_input < tmp_misobj_high;
    misobj_input_sel = misobj_input(index_sel_2);
    
    pd_2 = fitdist(misobj_input_sel,'normal');
    
    if abs(tmp_mu - pd_2.mu) < 0.1
        break
    else
        tmp_mu = pd_2.mu;
    end
end

misobj_fit = pd_2.mu;
% figure; histfit(misobj_input_sel,21,'normal')



end
