src = '/Volumes/G_PSYPLAFOR_methlab$/EEGManyLabs/data/prepAdvanced/';

dest = '/Volumes/BA891B30-33B9-41A7-B7D3-F67985C82D0F/VogelMachizawa2004/data/prepAdvanced/';
mkdir(dest)
% get all matching folders
d = dir(fullfile(src, 'sub-Vog*'));
d = d([d.isdir]);   % keep only directories

for i = 1:numel(d)
    
    disp(i)

    s = fullfile(src, d(i).name);
    des = fullfile(dest, d(i).name);

    copyfile(s, des);
    
end

%%
src = '/Volumes/G_PSYPLAFOR_methlab$/EEGManyLabs/data/prepICA/';

dest = '/Volumes/BA891B30-33B9-41A7-B7D3-F67985C82D0F/VogelMachizawa2004/data/prepICA/';
mkdir(dest)
% get all matching folders
d = dir(fullfile(src, 'sub-Vog*'));
d = d([d.isdir]);   % keep only directories

for i = 1:numel(d)
    
    disp(i)

    s = fullfile(src, d(i).name);
    des = fullfile(dest, d(i).name);

    copyfile(s, des);
    
end