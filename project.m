
clear all;
clc;
close all;

pkg load signal


function plotsignal(y,fs,header,showPlot="on")
    t=linspace(0,length(y)/fs,length(y));
    figure('name',header,"visible",showPlot);
    subplot(3,1,1);
    plot(t,y,'color',"b");
    xlabel("time domain");
    ylabel("Amp");
    
    mag=abs(fftshift(fft(y)));
    f=linspace(-fs/2,fs/2,length(mag));
    subplot(3,1,2);
    plot(f,mag,'color',"r");
    xlabel("Freq domain");
    ylabel("Magnitude");
    
    yangle=angle(fftshift(fft(y)));
    subplot(3,1,3);
    plot(f,yangle,'color',"g");
    xlabel("Freq domain");
    ylabel("Phase");
endfunction

function yb=bpf(ym,order,bw,fc,fs)
   [b,a]=butter(order,[fc-bw/2,fc+bw/2]/(fs/2));
   yb=filter(b,a,ym)*2;
endfunction

function ym = modulatesignal(y,fs,fc)
    t=linspace(0,length(y)/fs,length(y));
    ym=(y)'.*cos(2*pi*fc*t);
endfunction

function yd=demodulate(yb,fs,fc,dfc=0,dphi=0)
t=linspace(0,length(yb)/fs,length(yb));
yd=(yb)'.*cos(2*pi*(fc+dfc)*t+dphi);
endfunction

function yf=lpf(yd,order,BW,fs)
fc=BW/2;
[b,a]=butter(order,fc/(fs/2));
yf=filter(b,a,yd)*2;
endfunction



[y2,fs1]=audioread("fatha.wav");
[y1,fs2]=audioread("khaf.wav");
[y3,fs3]=audioread("fosalat.wav");


[y1,y2,y3]=deal(resample( y1(:,1),3,1),resample( y2(:,1),3,1),resample( y3(:,1),3,1));
d = max([length(y1) length(y2) length(y3)]);

padding=@(y,d) [y;zeros(d-length(y),1)];
[y1,y2,y3]=deal(padding(y1,d),padding(y2,d),padding(y3,d));
y=[y1,y2,y3];
fs=3*min(fs1,fs2,fs3);
%for yi=y, plotsignal(yi,fs,'Baseband signal'),endfor

plotsignal(y1,fs,'Baseband signal of first signal ');
plotsignal(y2,fs,'Baseband signal of second signal');
plotsignal(y2,fs,'Baseband signal of third signal');

fc1=10000;
fc2=30000;
fc3=50000;

ym1=modulatesignal(y1,fs,fc=10000);
ym2=modulatesignal(y2,fs,fc=30000);
ym3=modulatesignal(y3,fs,fc=50000);
ym=ym1+ym2+ym3;

plotsignal(ym,fs,'modulated signal');

signalsafterpassingbpf=[];
for fc=[fc1,fc2,fc3], signalsafterpassingbpf=[signalsafterpassingbpf,(bpf(ym,10,10000,fc,fs))'];,endfor
for ybx=signalsafterpassingbpf,plotsignal(ybx,fs,'bpf'),endfor





demodsignals=[];
fc=[fc1,fc2,fc3];
for i=1:length(fc), demodsignals=[demodsignals,(demodulate(signalsafterpassingbpf(:,i),fs,fc(i)))'];,endfor
for ydx=demodsignals,plotsignal(ydx,fs,'Demodulate without LPF'),endfor


yf=[];
for ydx=demodsignals, yf=[yf,(lpf(ydx',10,10000,fs))'];,endfor
for yfx=yf,plotsignal(yfx,fs,' LPF '),endfor

sound(yf(:,2),fs);

sound(yf(:,1),fs);
sound(yf(:,3),fs);

yd2=[];
phi=[deg2rad(10),deg2rad(30),deg2rad(90)];
for i=1:length(phi), yd2=[yd2,(demodulate(signalsafterpassingbpf(:,2),fs,fc(2),dfc=0,dphi=phi(i)))'];,endfor

yf2=[];
for ydx=yd2, yf2=[yf2,(lpf(ydx',10,10000,fs))'];,endfor

plotsignal(yf2(:,1),fs,'phi');
plotsignal(yf2(:,2),fs,'phi');
plotsignal(yf2(:,3),fs,'phi');

yd2=[];
df=[2,10];
for i=1:length(df), yd2=[yd2,(demodulate(signalsafterpassingbpf(:,2),fs,fc(2),dfc=df(i)))'];,endfor

yf2=[];
for ydx=yd2, yf2=[yf2,(lpf(ydx',10,10000,fs))'];,endfor
plotsignal(yf2(:,1),fs,'f');
plotsignal(yf2(:,2),fs,'f');


%sound(yf2(:,1),fs);
%sound(yf2(:,2),fs);















































  
