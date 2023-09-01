clear all; close all;
tstart=tic
% PLOT PDF
datalength=600;% length of one PSD line (sec) (600.30.240) (1800.60.480)
aa=30;% shifting every aa s
bb=240;% whole window bb s
cc=(datalength-bb)/aa+1;
year= num2str([2006:2019].','%04d');% year 
mmmm= num2str([01:12].','%02d');% month 
dddd= num2str([01:31].','%02d');% day 
kkkk= num2str([00:23].','%02d');% hour 
hours= 24;% change with counts of kkkk 
BHS= ['1';'2';'3';'4';'5';'6';'7'];% from 1 to 7
ENZ= ['E';'N';'Z'];% E or N or Z
%For data from 2006.10 to 2007.12, change fn: '-SH' into '-GH'.
totalpsd= 0;
%mat= '200608BHS1data'

%---------------------------main program----------------------------------%
% NLNM
[NLNM_data]=load('NLNM.txt');
period2=NLNM_data(:,1);
% convert to velocity
%NLNM=NLNM_data(:,2)+NLNM_data(:,3).*log10(period2)+20*log10(period2/2/pi);
%acc 
NLNM=NLNM_data(:,2)+NLNM_data(:,3).*log10(period2);
% NHNM
[NHNM_data]=load('NHNM.txt');
period1=NHNM_data(:,1);
% convert to velocity
%NHNM=NHNM_data(:,2)+NHNM_data(:,3).*log10(period1)+20*log10(period1/2/pi);
%acc 
NHNM=NHNM_data(:,2)+NHNM_data(:,3).*log10(period1);

for nn=1 %BHS
for ee=3 %ENZ
for yy=14 % which year (yy+2005)
for mm=8 % which month 
    totalpsd=0;
    %tictime1=tic
for dd=1:31%1:31 % which day 
for kk=10:21%1:24 % whole day hours 

        %fn1=['fish_TCDP_SAC/',year,'/',mm,'_',year,'/',mm,dddd(1,:),'_',year,'/T',year,mm,dddd(1,:),'-',kkkk(1,:),'00-BHS',BHS,'-SH',ENZ,'.sac']
        data_path1=['/DARRAY5/databases_sac/fish_TCDP_SAC'];
        data_path2=[year(yy,:),'/',mmmm(mm,:),'_',year(yy,:),'/',mmmm(mm,:),dddd(dd,:),'_',year(yy,:)];
        fn=['T',year(yy,:),mmmm(mm,:),dddd(dd,:),'-',kkkk(kk,:),'00-BHS',BHS(nn,:),'-SH',ENZ(ee,:),'.sac']; % read SAC filename 
        whole_path=[data_path1,'/',data_path2,'/',fn];
        %whole_path=['C:\Users\fish\Desktop\PSD\sac2\',fn];
        %tt=[year(yy,:),mmmm(mm,:),'-BHS',BHS(nn,:),'-',ENZ(ee,:)]; % title of the figure
        tt=[year(yy,:),mmmm(mm,:),'-BHS',BHS(nn,:),'-',ENZ(ee,:),'night']; % title of the figure
        fnpsd=['PSD-',tt]; %filename of figure(1)
        fnpsdr=['Resampled PSD-',tt]; %filename of figure(3)
        fnppsd=['PPSD-',tt]; %filename of figure(2)
        %fn='./20190729003807068.4.EHE.HL_rm.sac';
        %end

        %if exist(fn,'file')
        %[prey, header ]= sac2mat(fn);        
        if (exist(whole_path,'file') ~= 0)
           %load(whole_path); 
           [prey, header ]= sac2mat(whole_path);
           dt=header(1);
           npts=header(80);
           b=header(6);
           zero=0;
           for i=1:npts-4 %5 zeros then abandon
               if (prey(i) == 0 && prey(i) == prey(i+1) && prey(i+1) == prey(i+2) && prey(i+2) == prey(i+3) && prey(i+3) == prey(i+4))
                   zero=zero+1;
               end
           end
           if (zero == 0)

        %else
        %   disp('File not exist:'), disp(fn)
        %end
        
        for i=1:npts-1
            yd(1,i)=(prey(i+1)-prey(i))/dt;% differential to acceleration �L�����[�t��
        end
        %yd(1,npts)=yd(1,i);
        gg=int32(npts/(datalength/dt));
        
        for qq=1:gg  %(day)2*6:8*6 %(night)16*6:21*6
            if (qq==gg)
               y=yd(int32((qq-1)*datalength/dt):int32(qq*datalength/dt-1));
            else
               y=yd(int32((qq-1)*datalength/dt+1):int32(qq*datalength/dt)); 
            end
            %t= b : dt : (npts-1).*dt;
            % shifting every aa s
            shift=int32(aa/dt);
            % whole window bb s
            whole=int32(bb/dt);

            for i = 1:cc
                if (i==cc)
                amp=y((i-1)*shift:(i-1)*shift+whole-1);
                amp = amp *0.01;
                amp = amp-mean(amp);
                amp = detrend(amp);
                tp = tukeywin(whole,0.1)';
                amp = amp .* tp;
                [ psdx ,freq ] = sub_psd(amp,dt,length(amp));
                %psdall(i,25001) = NaN;
                %psdall(i,1:25000) = psdx;
                psdall(i,:) = psdx;
                else
                amp=y((i-1)*shift+1:(i-1)*shift+whole);
                % convert cm -> m
                amp = amp *0.01;
                % rmean
                amp = amp-mean(amp);
                % detrend
                amp = detrend(amp);
                % taper 10%
                tp = tukeywin(whole,0.1)';
                amp = amp .* tp;
                [ psdx ,freq ] = sub_psd(amp,dt,length(amp));
                psdall(i,:) = psdx;
                end
            end

            psdave(totalpsd+1,:) = mean(psdall);%averaged psd
    
    %------------------------resampling--------------------------------
            fff = log10(freq);
            ys = psdave(totalpsd+1,:);
            pt = 0.002;
            %pt = 0.1;
            fnew(totalpsd+1,:) = fff(end):pt:fff(2);
            sznew = length(fnew(totalpsd+1,:));
        
            for i = 1 : sznew
                item = fff < (fnew(totalpsd+1,i) + pt) & fff > (fnew(totalpsd+1,i) - pt);
                ysnew(totalpsd+1,i) = mean(ys(item));
            end
    
            %psdy((dd-1)*hours+kk,:) = ysnew;
            %psdp((dd-1)*hours+kk,:) = fnew;
            psdy(totalpsd+1,:) = ysnew(totalpsd+1,:);
            psdp(totalpsd+1,:) = fnew(totalpsd+1,:);    
            totalpsd=totalpsd+1;
            [ txt ] = file_psd(totalpsd,dd,kk,qq,ysnew); %built .mat file to save data 
            txttxt(totalpsd,:)=txt(1,:);
    
            %ii = (dd-1)*hours+kk;
        end
           else
               disp('File BROKEN:'), disp(fn)
           end
           else
           disp('File NOT EXIST:'), disp(fn)
        end
end
end
totalpsd

if exist('totalpsd','var') ~= 0 && (totalpsd ~= 0)
%    return;
%end


if exist('txttxt','var') ~= 0
sv=[fnpsd,'.mat'];
save(sv,'txttxt');
end

%toctime1=toc(tictime1)

%figure (1); % PSD
%for i=1:totalpsd
%    semilogx(freq,psdave(i,:));
%    hold on;
%end
%plot(period2,NLNM,'k')
%plot(period1,NHNM,'k')
%text(0.05,-91,'NHNM')
%text(0.05,-167,'NLNM')
%grid on
%title(['PSD-' tt])
%xlabel('Period (sec)')
%ylabel('PSD (dB)')
%ylim([-200 -80])
%xlim([0.01 1])
%print(fnpsd,'-dtiffn')
%close(1)
%tictime2=tic
figure(3) % PSD resampling
%semilogx(freq,psdave,'k-','linewidth',1)
%hold on;
for i=1:totalpsd
    str = num2str(i,'%01d');
    semilogx(10.^fnew(i,1:1002),ysnew(i,1:1002),'r-');
    text(10.^fnew(i,1000),ysnew(i,1000),str)
    hold on
end
grid on;
%title(['PSD resampling-' tt])
plot(period2,NLNM,'Color',[0.5 0.5 0.5])
plot(period1,NHNM,'Color',[0.5 0.5 0.5])
text(0.05,-167,'NLNM')
text(0.05,-91,'NHNM')
title(['PSD-' tt])
ylim([-200 -80])
xlim([0.01 1])
xlabel('Period (sec)')
ylabel('PSD (dB)')
%multitaper
%loglog(f,ym,'b-')
%resampling

print(fnpsd,'-dpdf','-r300')
print(fnpsd,'-djpeg','-r300')
close(3)
%toctime2=toc(tictime2)

%return
%----------------------------PPSD--------------------------------
svv=[fnppsd,'-tmp.mat'];
save (svv)  % Just a tmp file to make the program run faster.
clearvars -except svv; 
load (svv)
%tictime3=tic
[pdflinenum pdfcolnum]=size(psdy)

    for j=1:pdfcolnum
        k=1;
        dB_start=-80; %-40
        for dB=dB_start:-1:-200 
            pdfnum(k,j) = length(find((psdy(:,j) <= dB) & (psdy(:,j) > dB-1))); 
            pdffin(k,j) = pdfnum(k,j)/pdflinenum; % to calculate the percentage of PSD inside a 1dB bin.
            %figure (1);
            dBBB(k,j) = dB;
            psdppp(k,j) = psdp(1,j);
            k=k+1;
        end
        [max_num(j),max_idx(j)] = max(pdfnum(:,j));
        if  max_num(j) ~= 0
            MODE(j) = (max_idx(j))*(-1)+(dB_start+1);
        else
            MODE(j) = NaN;
        end
    end
%toctime3=toc(tictime3)
save (svv)  % Just a tmp file to make the program run faster.
clearvars -except svv;
load (svv)
%tictime4=tic
figure (2);
    for i=1:(200+dB_start+1)
        scatter(10.^psdppp(i,1:1002),dBBB(i,1:1002),20,pdffin(i,1:1002),'s','filled') %only use NO.1-NO.1002 point (because of resampling)
        hold on;
    end
set(gca,'xscale','log')
%colormap(jet)
cptcmap('GMT_rainbow');
caxis([0 0.3])
colorbar
ylabel(colorbar,'Probability')
hold on;
grid on;
plot(period2,NLNM,'Color',[0.5 0.5 0.5])
plot(period1,NHNM,'Color',[0.5 0.5 0.5])
plot(10.^psdppp(1,1:1002),MODE(1:1002),'k','LineWidth',0.5)
text(0.05,-167,'NLNM')
text(0.05,-91,'NHNM')
title(['PPSD-' tt])
ylim([-202 dB_start])
xlim([0.0095 1])
xlabel('Period (sec)')
ylabel('PSD (dB)')
%print(fnppsd,'-dtiffn')
print(fnppsd,'-dpdf','-r300')
print(fnppsd,'-djpeg','-r300')
close(2)
%toctime4=toc(tictime4)
delete (svv)
clearvars psdy psdp pdflinenum pdfcolnum pdfnum pdffin dBBB psdppp max_num max_idx MODE txttxt;
end
end
end
end
end
tend=toc(tstart)