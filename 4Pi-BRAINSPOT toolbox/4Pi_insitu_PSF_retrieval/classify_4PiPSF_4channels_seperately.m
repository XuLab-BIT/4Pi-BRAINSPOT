
function [ims_Zcal_ave_planes,index_record_Zplanes] = classify_4PiPSF_4channels_seperately(subregion_chs,ref_planes, empupil) 


%% Calculate the similarity between reference images and single molecules 

disp('Calculate the similarity between reference images and single molecules');
img_planes = double(subregion_chs);    %normalization image

num_img = size(img_planes,3);
num_ref = size(ref_planes,3);
num_ch = size(ref_planes,4);
imsz = empupil.imsz;  


similarity_in_ch1 = zeros(num_ref, num_img); 
similarity_in_ch2 = zeros(num_ref, num_img);
similarity_in_ch3 = zeros(num_ref, num_img);  
similarity_in_ch4 = zeros(num_ref, num_img); 


index_similarity = zeros(num_ref, num_img); %index


shift_row_ch1 = zeros(num_ref, num_img);
shift_col_ch1 = zeros(num_ref, num_img);
shift_row_ch2 = zeros(num_ref, num_img);
shift_col_ch2 = zeros(num_ref, num_img);
shift_row_ch3 = zeros(num_ref, num_img);
shift_col_ch3 = zeros(num_ref, num_img);
shift_row_ch4 = zeros(num_ref, num_img);
shift_col_ch4 = zeros(num_ref, num_img);

%parfor, XY shift seperately
parfor ii = 1 : num_img

    for jj = 1 : num_ref
             
        [shift1, shift2, tmpval1] = registration_in_each_channel(ref_planes(:,:,jj,1), img_planes(:,:,ii,1));
             
        if (abs(shift1) > 6 || abs(shift2) > 6 || tmpval1 < 0.4)
            continue;
        end
        
        %get X, Y shift in plane1
        shift_row_ch1(jj, ii) = shift1;
        shift_col_ch1(jj, ii) = shift2;
    
        [shift3, shift4, tmpval2] = registration_in_each_channel(ref_planes(:,:,jj,2), img_planes(:,:,ii,2));

        if (abs(shift3) > 6 || abs(shift4) > 6 || tmpval2 < 0.4)  
            continue;
        end
        
        %get X, Y shift in plane2
        shift_row_ch2(jj, ii) = shift3;
        shift_col_ch2(jj, ii) = shift4;
        
        %distance of two plane shift
        dist_shift_2planes = sqrt((shift1-shift3)^2 + (shift2-shift4)^2);
        if dist_shift_2planes > 2
            continue;
        end
        
        [shift5, shift6, tmpval3] = registration_in_each_channel(ref_planes(:,:,jj,3), img_planes(:,:,ii,3));

        if (abs(shift5) > 6 || abs(shift6) > 6 || tmpval3 < 0.4)  
            continue;
        end
        
        %get X, Y shift in plane3
        shift_row_ch3(jj, ii) = shift5;
        shift_col_ch3(jj, ii) = shift6;
        
        %distance of two plane shift
        dist_shift_2planes = sqrt((shift3-shift5)^2 + (shift4-shift6)^2);
        if dist_shift_2planes > 2
            continue;
        end
        
        [shift7, shift8, tmpval4] = registration_in_each_channel(ref_planes(:,:,jj,4), img_planes(:,:,ii,4));
        
        if (abs(shift7) > 6 || abs(shift8) > 6 || tmpval4 < 0.4)
            continue;
        end
        
        %get X, Y shift in plane4
        shift_row_ch4(jj, ii) = shift7;
        shift_col_ch4(jj, ii) = shift8;
        
        %distance of two plane shift
        dist_shift_2planes = sqrt((shift5-shift7)^2 + (shift6-shift8)^2);
        if dist_shift_2planes > 2
            continue;
        end        
        
        
        if (tmpval1+tmpval2+tmpval3+tmpval4)/4 < empupil.min_similarity
            continue;
        end
       
        
        similarity_in_ch1(jj, ii) = tmpval1;
        similarity_in_ch2(jj, ii) = tmpval2;
        similarity_in_ch3(jj, ii) = tmpval3;
        similarity_in_ch4(jj, ii) = tmpval4;
         
            
    end
   
end


%%
similarity_in_chs = (similarity_in_ch1 + similarity_in_ch2 + similarity_in_ch3 + similarity_in_ch4)./4; %similarity in 4Pi
for ii = 1 : num_img
    % Sort the similarity
    [sort_similarity, index_sort] =  sort(similarity_in_chs(:,ii),'descend'); 
    if (sort_similarity(1) == 0)
        continue;
    end
    % Determine single molecule image belong which reference image
    index_similarity(index_sort(1), ii) = 1;
    for jj = 2 : num_ref
        if (sort_similarity(jj) >= sort_similarity(1)-0.0) && (abs(index_sort(jj)-index_sort(1)) == 1)   
            index_similarity(index_sort(jj),ii) = 1;
        else
            break;
        end
    end
    
end


%% Updata average images 
disp('Update average images in 4 planes');

ims_Zcal_ave_planes = zeros(imsz,imsz,num_ref,num_ch);
index_record_Zplanes = zeros(num_ref,1);
for ii = 1 : num_ref
    index_selection = find(index_similarity(ii,:) == 1);
    sz_index = size(index_selection,2);
    if sz_index > empupil.bin_lowerBound   
        ims_planes_shift = zeros(imsz,imsz,sz_index,num_ch);
        for jj = 1 : sz_index

            ims_planes_shift(:,:,jj,1) = FourierShift2D(similarity_in_ch1(ii, index_selection(jj)) .* img_planes(:,:,index_selection(jj),1), [shift_row_ch1(ii, index_selection(jj)) shift_col_ch1(ii, index_selection(jj))]);
            ims_planes_shift(:,:,jj,2) = FourierShift2D(similarity_in_ch2(ii, index_selection(jj)) .* img_planes(:,:,index_selection(jj),2), [shift_row_ch2(ii, index_selection(jj)) shift_col_ch2(ii, index_selection(jj))]);
            ims_planes_shift(:,:,jj,3) = FourierShift2D(similarity_in_ch3(ii, index_selection(jj)) .* img_planes(:,:,index_selection(jj),3), [shift_row_ch3(ii, index_selection(jj)) shift_col_ch3(ii, index_selection(jj))]);
            ims_planes_shift(:,:,jj,4) = FourierShift2D(similarity_in_ch4(ii, index_selection(jj)) .* img_planes(:,:,index_selection(jj),4), [shift_row_ch4(ii, index_selection(jj)) shift_col_ch4(ii, index_selection(jj))]);

        end
        
        % average the images
        index_record_Zplanes(ii) = 1;
        for kk = 1: num_ch
            ims_Zcal_ave_planes(:,:,ii,kk) = mean(ims_planes_shift(:,:,:,kk),3);
        end
    end
end


%% remove the PSF too out-of-focus
Z_pos_ref = empupil.Z_pos + empupil.zshift;
for ii = 1 : num_ref
    Z_pos_sel = Z_pos_ref(ii);
    if Z_pos_sel < -0.6 || Z_pos_sel > 0.6 
        index_record_Zplanes(ii) = 0;
    end
end


%% Remove too far away Reassembled PSF

tmp_record_keep = find(index_record_Zplanes == 1);
size_tmp = length(tmp_record_keep);

for ii = ceil(size_tmp/2) : -1 : 2
    if tmp_record_keep(ii) > tmp_record_keep(ii-1)+3
        index_record_Zplanes(tmp_record_keep(ii-1)) = 0;
        tmp_record_keep(ii-1) = tmp_record_keep(ii);
    end
end

for ii = ceil(size_tmp/2) : size_tmp-1
    if tmp_record_keep(ii) < tmp_record_keep(ii+1)-3
        index_record_Zplanes(tmp_record_keep(ii+1)) = 0;
        tmp_record_keep(ii+1) = tmp_record_keep(ii);
    end
end 
