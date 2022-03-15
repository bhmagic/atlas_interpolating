clear

csv_name = [pwd '/atlas_info_KimRef_FPbasedLabel_v2.7.csv'];

csv = readmatrix(csv_name,'OutputType','string');
csv_d = readmatrix(csv_name,'OutputType','double');


index_id = 1;
index_parent_id = 8;
index_name = 2;
index_acronym = 3;
index_structure_order = 7;



T = readtable(csv_name);

id = table2array(T(:,index_id));
parent_id = table2array(T(:,index_parent_id));
name = table2array(T(:,index_name));
acronym = table2array(T(:, index_acronym));
structure_order = table2array(T(:, index_structure_order));


bad_parent_idx = ~ismember(parent_id,id);
bad_parent = parent_id(bad_parent_idx)
bad_parent_kid = id(bad_parent_idx)

idx = find(~isnan(id));
[~,p_idx]=ismember(parent_id,id);


% [c, ia, ic] = unique(name);
% c(find(accumarray(ic,1)>2))
% [c, ia, ic] = unique(acronym);
% c(find(accumarray(ic,1)>2))

% G = digraph(p_idx(3:end), idx(3:end), 1, name);

depth_level = ones(size(name));
new_d = depth_level;
new_d(idx(3:end)) = depth_level(p_idx(3:end))+1;

while ~all(new_d==depth_level)
    depth_level = new_d;
    new_d(idx(3:end)) = depth_level(p_idx(3:end))+1;
    new_d = max([depth_level, new_d],[],2);
    
end



xxx = [depth_level(~bad_parent_idx), depth_level(p_idx(~bad_parent_idx)) ]-1;
yyy = [idx(~bad_parent_idx), idx(p_idx(~bad_parent_idx)) ];

figure;scatter(depth_level-1,idx);hold;
text(depth_level+0.1-1,idx,name);
plot(xxx',yyy');

ylim([0,64]);
xlim([0,24]);
