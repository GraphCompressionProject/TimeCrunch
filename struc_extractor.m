input_dir = 'path/to/data/directory/';
ext = 'txt'
flist = dir([input_dir,strcat('*.',ext)])
numfiles = numel(flist)
for i = 1:numfiles
	fn = [input_dir,flist(i).name];
	run_structureDiscovery(fn);
end	
