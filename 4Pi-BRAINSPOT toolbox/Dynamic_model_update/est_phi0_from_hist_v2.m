%% estimate the cavity phase parameter from measured histgraph

function [phi0_guess_keep,phi0_fit,totmask] = est_phi0_from_hist_v2(phi0_guess,z_guess,mcc_val)


%% rejection
ccmask = mcc_val < 0.6; 
zmask = z_guess > 0.5 | z_guess < -0.5; 
totmask = ccmask | zmask;
phi0_guess_keep = phi0_guess(~totmask);

%% parameter estimation

mean_phi0 = angle(mean(exp(1i.*phi0_guess_keep)));
phi0_input = [phi0_guess_keep-2*pi; phi0_guess_keep; phi0_guess_keep+2*pi];

range_cal = pi;
tmp_phi0_low =  mean_phi0 - range_cal;
tmp_phi0_high = mean_phi0 + range_cal;

index_sel = phi0_input >= tmp_phi0_low & phi0_input <= tmp_phi0_high;
phi0_input_sel = phi0_input(index_sel);

pd = fitdist(phi0_input_sel,'normal');
tmp_mu = pd.mu;

for ii = 1 : 5

    tmp_phi0_low =  tmp_mu - range_cal;
    tmp_phi0_high = tmp_mu + range_cal;
    
    index_sel_2 = phi0_input > tmp_phi0_low & phi0_input < tmp_phi0_high;
    phi0_input_sel = phi0_input(index_sel_2);
    
    pd_2 = fitdist(phi0_input_sel,'normal');
    
    if abs(tmp_mu - pd_2.mu) < 0.1
        break
    else
        tmp_mu = pd_2.mu;
    end
end

phi0_fit = pd_2.mu;


end
