%A helper function to read a GPS log text file
%Tries to handle both , and . as the decimal separator
%@params fileIn full path to the file to be read
%@params separator the column separator, defaults to tab '\t'
%@params headerLines, number of lines at the top of the file to consider header
%@return struct data with fields header and data
function data = readLog(fileIn,separator,headerLines)
	if ~exist('separator','var')
		separator = '\t';
	end
	if ~exist('headerLines','var')
		headerLines = 1;
	end
	data = struct();
	%Read header lines
	tempHeader = struct();
	fh = fopen(fileIn,'r');
	for i = 1:headerLines
		tempHeader(i).line = strsplit(fgetl(fh),separator);
	end
	data.header = tempHeader;
	
	try
		%Get a line of data with either . or , in it to check which is used as a decimal separator
		temp = strsplit(fgetl(fh),separator);
		while sum(cellfun(@(x) ~isempty(findstr(x,'.')),temp)) == 0 && sum(cellfun(@(x) ~isempty(findstr(x,',')),temp)) == 0
			temp = strsplit(fgetl(fh),separator);
		end
		fclose(fh);
	
		if sum(cellfun(@(x) ~isempty(findstr(x,'.')),temp)) > 0
			%Likely . as the decimal separator
			data.data = dlmread(fileIn,separator,headerLines,0);
		else
			%Likely , as the decimal separator, read all bytes into memory, replace all commas with dots
			fh = fopen(fileIn,'r');
			fileBytes = uint8(fread(fh));	%Read data into memory
			fclose(fh);
%             keyboard;
%             if exist('OCTAVE_VERSION', 'builtin')
%                 fileBytes(fileBytes == 0x2c) = uint8(0x2e);	%Replace all commas (0x2c) with dots (0x2e)
%             else
                fileBytes(fileBytes == uint8(hex2dec('2c'))) = uint8(hex2dec('2e'));	%Replace all commas (0x2c) with dots (0x2e)
%             end
			%Write data into a temporary file
			tmpName = [tempname '.txt'];	%Get temporary file name
			fh = fopen(tmpName,'w');
			fwrite(fh,fileBytes);	%write bytes to the temporary file
			fclose(fh);
			data = readLog(tmpName,separator,headerLines);	%Read the file with this function
			delete(tmpName);	%Remove the temporary file
		end	
	catch
	end
	
