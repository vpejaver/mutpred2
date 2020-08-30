% filter_twosided.m

% data is filtered by columns

function fdata=filter_twosided(data,filt_onesided)

one=ones(size(data));
scale1=filter(filt_onesided,1,one);
scale2=flipud(filter(filt_onesided,1,flipud(one)));
scale=(scale1+scale2)/2;

fdata1=filter(filt_onesided,1,data);
fdata2=flipud(filter(filt_onesided,1,flipud(data)));
fdata=(fdata1+fdata2)/2;
fdata=fdata./scale;

return;
