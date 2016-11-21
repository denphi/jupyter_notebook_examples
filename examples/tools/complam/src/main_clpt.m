function main_clpt (infile)

try
   lib = rpLib(infile);
   clpt(lib);
catch
   lastmsg = lasterr
   fid = fopen('errorDetected.txt','wt');
   fprintf(fid,'1\n');
   fprintf(fid,'%s\n',lastmsg);
   fclose(fid);
   rpLibPutString(lib,'output.log',lastmsg,1);
   rpLibResult(lib,1);
end
