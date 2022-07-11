function [All_net, L_net, R_net, homotopic, interhem] = ...
    net_mask_ind_reduced(NET, i_mask)

n = .5 + sqrt(.25 + 2*length(i_mask));    % # nodes n(n-1)/2 = l

% Find Within and Between Network mask
n_NET = length(NET.Index_reduction);  % # networks

All_net.mask = cell(n_NET, n_NET);

mask0 = zeros(n, n);
for in = 1 : n_NET
    for iin = in : n_NET
        mask_net = mask0;
        mask_net(NET.Index_reduction{in}, NET.Index_reduction{iin}) = 1;

        mask_net = triu(mask_net, 1);
        i_net    = find(mask_net);

        [~, All_net.mask{in, iin}] = intersect(i_mask, i_net);
    end
end

% Find only Left hemisphere network mask

if isstring(NET.LEFT_RIGHT)
    ind_L = find(strcmp(NET.LEFT_RIGHT, 'L'));
else
    ind_L = find(NET.LEFT_RIGHT == 'L');
end
L_net.mask = All_net.mask(ind_L, ind_L);
L_net.lab  = NET.Names(ind_L);

% Find only Right hemisphere network mask

if isstring(NET.LEFT_RIGHT)
    ind_R = find(strcmp(NET.LEFT_RIGHT, 'R'));
else
    ind_R = find(NET.LEFT_RIGHT == 'R');
end
R_net.mask = All_net.mask(ind_R, ind_R);
R_net.lab  = NET.Names(ind_R);    

% Find Homotopic networks mask

homotopic.mask  = cell(length(ind_L), 1);
homotopic.lab   = cell(length(ind_L), 1);
for in = 1 : length(ind_L)
    il = ind_L(in);
    ir = il + 1;
    
    mask_net = mask0;
    mask_net(NET.Index_reduction{il}, NET.Index_reduction{ir}) = 1;
    
    mask_net = triu(mask_net, 1);
    i_net    = find(mask_net);
    
    [~, homotopic.mask{in}] = intersect(i_mask, i_net);
    
    name = char(NET.Names(il));
    if il < 18
        homotopic.lab{in} = name(2 : end);
    else
        homotopic.lab{in} = name(1 : end-5);
    end
end
homotopic.lab = categorical(homotopic.lab);

% Find Interhemispheric (non-homotopic)  mask
if isstring(NET.LEFT_RIGHT)
    ind_none      = find(strcmp(NET.LEFT_RIGHT, 'none'));
else
    ind_none      = find(NET.LEFT_RIGHT == 'none');
end
ind_Ln        = union(ind_none, ind_L);
ind_Rn        = union(ind_none, ind_R);
interhem.mask = All_net.mask(ind_Ln, ind_Rn);
interhem.lab  = cell(length(ind_Ln), length(ind_Rn));

for in = 1 : length(ind_Ln)
    for iin = 1 : length(ind_Rn)
        if in == iin
            interhem.mask{in, iin} = [];
        elseif in > iin
            interhem.mask{in, iin} = All_net.mask{ind_Rn(iin), ind_Ln(in)};
        end
        
        if in ~= iin
            name1 = char(NET.Names(ind_Ln(in)));
            if (in > 26) && (in ~= 31)
                name1 = name1(1 : end-5);
            end
            name2 = char(NET.Names(ind_Rn(iin)));
            if (iin > 26) && (iin ~= 31)
                name2 = name2(1 : end-5);
            end
            
            interhem.lab{in, iin} = [name1, '-', name2];
        end
    end
end

% Identification of each pair of networks
% H: homotopic
% L: left
% R: right
% n: none

ID = cell(n_NET);
ID(ind_L, ind_L) = {'L'};
ID(ind_R, ind_R) = {'R'};

for in = 1 : length(ind_L)
    il = ind_L(in);
    ir = il + 1;
    ID{il, ir} = 'H';
end

for in = 1 : n_NET
    for iin = 1 : n_NET
        if iin < in
            ID{in, iin} = '';
        elseif isempty(ID{in, iin})
            ID{in, iin} = 'none';
        end
    end
end

ID = categorical(ID);

All_net.ID = ID;