function dtt = loadDTT(file_name)

file.name = file_name;

[file.path file.base file.ext] = fileparts(file.name);

% file.mat = [file.path '/' file.base '.mat'];
file.mat = ['./' file.base '.mat'];
if ~exist(file.mat, 'file')
    dtt = DttData(file.name,'verbose');
    save(file.mat,'dtt');
else
    load(file.mat,'dtt');
end

dtt.channels

end
