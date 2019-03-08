function bool = fileExists(fileName)
% 
% FILEEXISTS Check if file exists.
% fileExists(filename)
% 
% Return true if the file exists.
% 
if ispc
    import java.io.*;
    a=File(fileName);
    bool=a.exists();
elseif isunix
    bool=exist(fileName,'file')==2;
end

end