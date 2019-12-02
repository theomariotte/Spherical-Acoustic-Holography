function out_file_name = getOutFileName(Path,fname,ext)

if ~exist(Path,'dir')
    warning('Path does not exist. Creating new directory.');
    mkdir(Path);
end

if ext(1) ~= '.'
   ext = ['.' ext]; 
end


fID = 1;
fname2save = [fname sprintf('_%02d',fID)];

while(exist([Path fname ext],'file'))
    fID = fID + 1;
    fname2save = [fname sprintf('_%02d',fID)];
end


out_file_name =[Path fname2save ext];

end