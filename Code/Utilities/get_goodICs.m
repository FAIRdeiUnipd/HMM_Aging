function [comp_network_names,slctIC_origIC] = get_goodICs()

comp_network_names = {
    'VIS',      [9 31 32 42 51 83 84 118 152 153 157];      % Visual network 
    'SMN',      [1 6 3 4];                                  % Sensory-motor 
    'AUD',      [11];                                       % Auditory network
    'CON',      [22 35 110];                                % Cingulo-opercularis network 
    'DAN',      [77 81 113];                                % Dorsal-attention network
    'FPN',      [52 105 109 126];                           % Fronto-parietal network 
    'DMN',      [33 45 102 39 117 15];                      % Default-mode network
    'CCN',      [46 55 56 61 68 72 82 98 107];              % Cognitive-Control Net
    'FRN',      [19 124 145];                               % Frontal nodes, possibly hubs
%     'LAN',      [122]                                       % Language network
    'CER',      [14]                                        % Cerebellum
    'BG',       [5]};                                       % Basal ganglia  

slctIC_origIC = [comp_network_names{:,2}];
