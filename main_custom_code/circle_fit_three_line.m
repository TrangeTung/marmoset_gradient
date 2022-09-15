function [RADIUS,CENTER_loc] = circle_fit_three_line(frame,pupil_side)
possible = [0.01 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 1];
circular_sum = zeros(size(pupil_side));
line = circular_sum;
for p = 1:numel(possible) 
    for side = 1:2
        circular_plus = pupil_side(:,:,side);
        [source_row,source_col] = find(ones(size(frame))==1);
        [reference_row,reference_col] = find(circular_plus==1);
        circular_sum(:,:,side) = circular_plus;
        reference_rep_sub = permute(repmat([reference_row,reference_col],[1 1 numel(source_row)]),[1 3 2]);
        source_rep_sub = permute(repmat([source_row,source_col],[1,1,numel(reference_row)]),[3 1 2]);
        dif = reference_rep_sub-source_rep_sub;
        distance = sqrt(dif(:,:,1).*dif(:,:,1) + dif(:,:,2).*dif(:,:,2));
        % distance variation
        variation = std(distance,0,1);
        variation_up = sort(variation);
        center_var_index = find(variation < variation_up(ceil(numel(variation)*possible(p))));
        center_var = zeros(size(frame));
        center_var(center_var_index) = 1;
        % distance mean
        mean_dist = mean(distance,1);
        line(:,:,side) = center_var;
        dis(side,:) = mean_dist;
    end
    circular_leftright = sum(circular_sum,3);
    dis_diff = abs(dis(1,:)-dis(2,:));
    dis_diff_up = sort(dis_diff);
    center_dis_diff_index = find(dis_diff < dis_diff_up(ceil(numel(dis_diff)*possible(p))));
    center_dis_diff = zeros(size(frame));
    center_dis_diff(center_dis_diff_index) = 1;
    
    center_possible = line(:,:,1).*line(:,:,2).*center_dis_diff;
    % possible center
    [center_possible_row,center_possible_col] = find(center_possible==1);
    center_possible_index = find(center_possible==1);
    if ~isempty(center_possible_index)
        break;
    end
end
% figure;imshow(line(:,:,1)*1+line(:,:,2)*4+center_dis_diff*8,[])
%% cicular fitting

distance_min = zeros(length(center_possible_row),1);
for iter = 1:length(center_possible_row)
    % background = fitting_pupil_index;
    radius = mean_dist(center_possible_index(iter));
    background = rgb2gray(insertShape(zeros(size(frame)),'circle',[center_possible_col(iter),...
        center_possible_row(iter) radius],'LineWidth',1,'Color','w'));
    % min_distance
    [source_row,source_col] = find(background~=0);
    [reference_row,reference_col] = find(circular_leftright==1);
    reference_rep_sub = permute(repmat([reference_row,reference_col],[1 1 numel(source_row)]),[1 3 2]);
    source_rep_sub = permute(repmat([source_row,source_col],[1,1,numel(reference_row)]),[3 1 2]);
    dif = reference_rep_sub-source_rep_sub;
    distance = sqrt(dif(:,:,1).*dif(:,:,1) + dif(:,:,2).*dif(:,:,2)).*repmat(background(background~=0)',[length(reference_row),1]);
    distance_min(iter,1) = mean(min(distance,[],2));
end
% center_loc & radius
center_index = find(distance_min == min(distance_min));
center_index(2:end) = [];
RADIUS = mean_dist(center_possible_index(center_index));
CENTER_loc = [center_possible_row(center_index),center_possible_col(center_index)];

end