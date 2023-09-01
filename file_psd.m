function [ txt ] = file_psd(totalpsd,dd,kk,qq,ysnew)

for i=1
    txt(i,1) = round(totalpsd);
    txt(i,2) = round(dd);
    txt(i,3) = round(kk-1);
    txt(i,4) = round(qq*10-10);
    txt(i,5:1006) = roundn(ysnew(totalpsd,1:1002),-3);
end