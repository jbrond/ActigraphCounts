%Helper function to remove . and ..
function fList = getFilesAndFolders(pathIn)
	fList = dir(pathIn);
	fNames = {fList(:).name};
	removeIndice = cellfun(@(x) strcmp(x,'.') | strcmp(x,'..'),fNames);
	fList(removeIndice) = [];
end