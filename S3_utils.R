# (C) Patrick Cahan 2012-2015
# utilitiy functions

s3_get_par<-function# get files in parallel
(bucket,
  path,
  fnames){
  tfname<-'tmpfile.txt';
  write.table(fnames,tfname, col.names=FALSE, row.names=FALSE, quote=FALSE);
  cmd<-paste("cat ", tfname, " | parallel aws s3 cp s3://", bucket, "/", path, "/{} ./", sep='');
  system(cmd);
}

utils_unpack<-function# gzip and untar
(fname){
  cmd<-paste("gzip -d ", fname, sep='');
  system(cmd);
  fname<-strsplit(fname, ".gz")[[1]][1];
  cmd<-paste("tar -xvf ", fname, sep='');
  system(cmd);
}

s3_get<-function
(dir,
 target, # name of target file
 bucket="pcahanrnaseq"
){
  fpath<-paste(bucket, "/",dir,sep='');
  cmd<-paste("aws s3 cp s3://", bucket, "/", dir, "/", target, " ./", sep='');
  system(cmd);  
}

s3_put<-function
(dir,
 target, # name of target file
 bucket="pcahanrnaseq"
){
  fpath<-paste(bucket, "/",dir,sep='');
  cmd<-paste("aws s3 cp ",target," s3://", bucket, "/", dir, "/", sep='');
  cat(cmd,"\n");
  system(cmd);  
}

utils_log<-function
#### make a logging string
(str, ### str to log
  prepend="" ### string to preprend
 ){
  lstring1<-paste(prepend, "[",format(Sys.time(), "%Y_%b_%d_%H:%M:%S"),"]\t", sep='');
  paste(lstring1, str, "\n", sep='');
  ### output str
}

s3_file_append<-function
#### add a line to a file on S3, always start with date
(dir, ### directory on bucket
 target, ### filename 
 str, ### str to log
 bucket="pcahanrnaseq"
 ){
  # fetch file from S3
  s3_get(dir, target, bucket);
  # add a line
  cat(utils_log(str), file=target, append=TRUE);
  # send file back
  s3_put(dir, target, bucket);
}

utils_tarcompress_dir<-function # tars and compresses (gzip) dir
(dir){
  fname<-paste(dir,".tar.gz", sep='');
  cmd<-paste("tar -cvzf ", fname," ",dir, sep='');
  system(cmd);
  fname;
  ### return name of tarred compressed file
}

utils_strip_fname<-function # reduces full path to filename
(str){
  a<-strsplit(str, "/")[[1]];
  a[length(a)];
}

utils_loadObject<-function
### loads an R object when you don't know the name
(fname
 ### file name
){
  x<-load(fname);
  get(x);
}

utils_count_comments<-function
(fname){
  cmd<-paste("grep \"#\" ", fname,sep='');
  x<-system(cmd, intern=T);
  length(x);
}

utils_stripwhite<-function
### strip whitespace from a string
(string
 #### string
 ){
  gsub("^\\s+|\\s+$", "", string)
}

utils_myDate<-function
### print date
()
{
  format(Sys.time(), "%b_%d_%Y");
}