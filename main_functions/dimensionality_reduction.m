%% ===================================================================== %%
%                        DIMENSIONALITY REDUCTION                         %
%  ===================================================================== %%
% 
% inputs: 
% - Coord: [N x 3] spatial coordinates of original data
% - NET_info: stucture with 3 fields:
%     1. Names: (Categorical [n_net x 1]) networks names
%     2. Index: (Cell [n_net x 1]) index of parcels in each network
%     3. LEFT_RIGHT (Categorical [n_net x 1]) L=left, R=right
% - sub_start: index of the first subcortical network
% 
% output:
% - NET_reduction: NET_info with additional fields:
%     1. Index_reduction: (Cell [n_net 1]) index of reduced parcels in each
%     network
%     2. final_index: (Cell [n_red, 1]) index of original parcellation
%     clusterized in each reduced parcel
%     3. net: [n_orig x 1] numerical vector indicating the network belongings
%     for each original parcel
%     4. net_reduced: [n_red x 1] numerical vector indicating the network belongings
%     for each reduced parcel
%     5. Names2: (Categorical [n_net_noLR, 1]) name of ecah network (without
%     L or R
%     6. Tick: [n_net_noLR, 1] numercial vector indicating the ticks for
%     original FC matrix plot
%     7. Tick_reduced: [n_net_noLR, 1] numercial vector indicating the ticks for
%     reduced FC matrix plot
%  ===================================================================== %%


function NET_reduction = dimensionality_reduction(Coord, NET_info, sub_start)

addpath('external/BCT')
addpath('utility')
addpath(genpath('external/simpleBrainSurface/'))

n_NET = length(NET_info.Index);

NET_info.Index_reduction = cell(n_NET, 1);

res_cl = cell(n_NET, 1);

N_orig = size(Coord, 1);
final_index = cell(N_orig, 1);


ik = 1;
n  = 0;
for in = 1 : n_NET
    idx_tmp   = NET_info.Index{in};
    
    if length(idx_tmp) <= 3
        n = n + length(idx_tmp);
        
        NET_info.Index_reduction{in} = (ik : (ik + length(idx_tmp) - 1))';
        for ii = 1 : length(idx_tmp)
            final_index{ik} = idx_tmp(ii);
            ik = ik + 1;
        end
    else
        coord_tmp = Coord(idx_tmp, :);
        Z  = linkage(coord_tmp, 'ward');
        
        res_cl{in} = evalclusters(coord_tmp, 'linkage', 'silhouette', ...
            'KList', [3 4 5]);
        
        [~, max_k] = max(res_cl{in}.CriterionValues);
        
        if max_k == 1
            Cl = cluster(Z,'Maxclust', 3); 
            NET_info.Index_reduction{in} = (ik : (ik + 2))';
            for ic = 1 : 3
                val = idx_tmp(find(Cl == ic));
                final_index{ik} = val;
                ik = ik + 1;
            end
            
            n = n + 3;
        elseif max_k == 2
            Cl = cluster(Z,'Maxclust', 4); 
            NET_info.Index_reduction{in} = (ik : (ik + 3))';
            for ic = 1 : 4
                val = idx_tmp(find(Cl == ic));
                final_index{ik} = val;
                ik = ik + 1;
            end
            
            n = n + 4;
        else
            Cl = cluster(Z,'Maxclust', 5); 
            NET_info.Index_reduction{in} = (ik : (ik + 4))';
            for ic = 1 : 5
                val = idx_tmp(find(Cl == ic));
                final_index{ik} = val;
                ik = ik + 1;
            end
            
            n = n + 5;
        end
    end
end

final_index = final_index(1 : ik-1);

NET_info.final_index = final_index;

Names     = NET_info.Names(1 : 19);
Names(19) = categorical("subcort");

Names2 = cell(10, 1);

net     = nan(343, 1);
net_red = nan(n, 1);

for in = 1 : 2 : sub_start
    if in < sub_start
        
        net([NET_info.Index{in}; NET_info.Index{in+1}]) = (in+1)/2;
        net_red([NET_info.Index_reduction{in}; NET_info.Index_reduction{in+1}]) ...
            = (in+1)/2;
    
        name_tmp = char(Names(in));
        Names2{(in+1) / 2} = name_tmp(2 : end);
    else
        
        net(NET_info.Index{in} : end) = (in+1)/2;
        net_red(NET_info.Index_reduction{in} : end) = (in+1)/2;
        
        Names2{end} = 'subcort';
    end
end

Names2 = categorical(Names2);

NET_info.net         = net;
NET_info.net_reduced = net_red;
NET_info.Names2      = Names2;

net_tick        = zeros(length(Names2), 1);
net_tick_red    = zeros(length(Names2), 1);

ik  = 0;
iik = 0;
for in = 1 : length(Names2)
    net_tick(in)     = ik + length(find(net == in))/2;
    net_tick_red(in) = iik + length(find(net_red == in))/2;
    
    ik  = ik + length(find(net == in));
    iik = iik + length(find(net_red == in));
end

NET_info.Tick         = net_tick;
NET_info.Tick_reduced = net_tick_red;


%% Plot Coordinates

% New Coordinates
Coord_reduced = nan(n, 3);

for in = 1 : n
    Coord_reduced (in, :) = mean(Coord(final_index{in}, :), 1);
end

Community_colors1 = my_color_map; %lines((max(Community_all)-2) / 2 + 2);
Community_colors1 = Community_colors1([1 : 3, 5 : end], :);

par.range = [.3 .7];

figure;
h = simpleBrainSurface(par, .5, 3);
hold on

for inet = 1 : max(net_red)
    
    idx = find(net_red == inet);
    
    plot3(Coord_reduced(idx, 1), Coord_reduced(idx, 2), ...
        Coord_reduced(idx, 3), 'o', 'MarkerFaceColor', Community_colors1(inet, :), ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 10)
end
    
hold off


%% Save Results

NET_reduction = NET_info; 

%% Remove folders from path

rmpath('external/BCT')
rmpath('utility')

rmpath(genpath('external/simpleBrainSurface/'))
