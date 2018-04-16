function writeToLog(file, text)
% write stuff to log file

fid = fopen(file, 'a');

fprintf(fid, '[%s]\n%s\n', datestr(now, 0),text);

fclose(fid);

end