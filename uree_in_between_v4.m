clear

original_anno_map = niftiread('uree_on_allen_20.nii');

upscaliing_factor = 5;

searching_range = 16;


original_anno_map = padarray(original_anno_map,[searching_range, searching_range, 0], 0);


radii_mask = inf(searching_range.*2+1);
[xxx,yyy] = ind2sub(size(radii_mask),find(radii_mask));

alt_ind = sub2ind([size(original_anno_map,1),size(original_anno_map,2)], xxx, yyy);
alt_ind = alt_ind - alt_ind(ceil(length(alt_ind)./2));

xxx = xxx - (searching_range+1);
yyy = yyy - (searching_range+1);
radii_mask = (xxx.*xxx + yyy.*yyy).^0.5;

new_map = {};

parfor ii = 1:size(original_anno_map,3)-1
% parfor ii =1:32
    og_up = original_anno_map(:,:,ii);
    og_down = original_anno_map(:,:,ii+1);

    missmatch_index = og_up == og_down;
    global_not_missing_locations = find(~missmatch_index);
    matching_distances_up = zeros(size(global_not_missing_locations));
    matching_distances_down = zeros(size(global_not_missing_locations));
    matching_index_up = zeros(size(global_not_missing_locations));
    matching_index_down = zeros(size(global_not_missing_locations));

    for jj = 1:length(global_not_missing_locations)

        local_masked_index = global_not_missing_locations(jj) + alt_ind;

        cd_1 = (og_up(global_not_missing_locations(jj)) == og_down(local_masked_index));
        cd_2 = og_down(global_not_missing_locations(jj)) ~= og_up(local_masked_index) & og_down(global_not_missing_locations(jj)) ~= og_down(local_masked_index);
        cd_3  = og_up(global_not_missing_locations(jj)) == og_up(local_masked_index) & og_up(local_masked_index) == og_down(local_masked_index);
        valid_mask = double(cd_1 | cd_3);
        if sum(valid_mask)==0
            valid_mask = double(cd_2);
        end

        [min_val, min_ind] = min(radii_mask./valid_mask);
        matching_distances_up(jj) = min_val;
        matching_index_up(jj) = local_masked_index(min_ind);


        cd_1 = (og_down(global_not_missing_locations(jj)) == og_up(local_masked_index));
        cd_2 = og_up(global_not_missing_locations(jj)) ~= og_up(local_masked_index) & og_up(global_not_missing_locations(jj)) ~= og_down(local_masked_index);
        cd_3 = og_down(global_not_missing_locations(jj)) == og_down(local_masked_index) & og_up(local_masked_index) == og_down(local_masked_index);

        valid_mask = double(cd_1 | cd_3);
        if sum(valid_mask)==0
            valid_mask = double(cd_2);
        end

        [min_val, min_ind] = min(radii_mask./valid_mask);
        matching_distances_down(jj) = min_val;
        matching_index_down(jj) = local_masked_index(min_ind);

    end


    z_offseting = (ii-1).*upscaliing_factor ;

    new_map{1,1,ii} = zeros([size(original_anno_map,1),  size(original_anno_map,2),  upscaliing_factor]);
    new_map{1,1,ii}(:,:,1) = og_up;
    temp = zeros(size(og_up));
    temp(missmatch_index) = og_up(missmatch_index);
    for jj = 2:upscaliing_factor
        temp2 = zeros(size(og_up));
        choice_flag = matching_distances_up.*(jj-1) < matching_distances_down.*(upscaliing_factor-jj+1);
        %         choice_flag = matching_distances_down.*(jj-1) > matching_distances_up.*(upscaliing_factor-jj+1);
        temp2(global_not_missing_locations(choice_flag)) = og_up(global_not_missing_locations(choice_flag));
        temp2(global_not_missing_locations(~choice_flag)) = og_down(global_not_missing_locations(~choice_flag));

        new_map{1,1,ii}(:,:,jj) = temp + temp2;
    end




end
new_map{1,1,size(original_anno_map,3)} = zeros([size(original_anno_map,1),  size(original_anno_map,2),  upscaliing_factor]);
new_map{1,1,size(original_anno_map,3)}(:,:,1) = original_anno_map(:,:,size(original_anno_map,3));



new_map = cell2mat(new_map);

new_map = new_map(searching_range+1:end-searching_range, searching_range+1:end-searching_range,:);



niftiwrite(new_map , 'uree_interpolate_v4.nii');


