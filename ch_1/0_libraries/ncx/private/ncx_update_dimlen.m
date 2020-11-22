function ncx_update_dimlen(theHandle,file,varname,dimName)

dim  = n_vardim(file,varname,dimName);
set(theHandle,'string',dim);

