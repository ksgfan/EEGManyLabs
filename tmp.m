

%%
src = '/Volumes/G_PSYPLAFOR_methlab$/EEGManyLabs/data/prepDirect'

dest = '/Volumes/BA891B30-33B9-41A7-B7D3-F67985C82D0F/VogelMachizawa2004/data/prepDirect/';
mkdir(dest)
% get all matching folders
d = dir(fullfile(src, '*'));
d(1:3) = [];

% for i = 112:217%numel(d)
for i = 112:114%numel(d)
    
    disp(d(i).name)

    s = fullfile(d(i).folder, d(i).name);
    
    des = fullfile(dest, d(i).name);
    
    mkdir(des)
    copyfile(s, des);

    
end
