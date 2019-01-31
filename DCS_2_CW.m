%% Preamble
% This code uses CW error signals to correct offset flucuations then
% resamples using DDFG. All data files must end in 'A' or 'B' e.g.
% 2018.07.26_17.25.04_1.1.A.bin

plots = 1*[1 1 1 1]; % 1 = create diagnostic plot, 0 = don't
t_delay = 0*10^-9; % Delay the correction signal by this amount in s
S = 4000*10^6; % sample rate in S/s
r = 20; % decimation factor
L = 10^7;
%% Load raw data
[filename, path] = uigetfile({'*.txt;*.bin' 'Text file or binary document';...
    '*.txt' 'Text file'; '*.bin' 'Binary document'},'Select BOTH raw data files',...
    'R:\DiodeDualComb\Data','MultiSelect','on');

% search for B file: if find B file then get name else use A name.
if size(filename,2) == 2
    filenameA = fullfile(path,filename{1});
    filenameB = fullfile(path,filename{2});
else
    filenameA = fullfile(path,filename);
    filenameB = filenameA;
end

if strcmp('.bin',filenameA(end-3:end)) == 1
    fidA = fopen(filenameA,'r');
    fidB = fopen(filenameB,'r');
    [ch1,~] = fread(fidA,[L,1],'uint16');
    [ch2,~] = fread(fidB,[L,1],'uint16');
    fclose(fidA);
    fclose(fidB);
    ch1 = 2.5*ch1/2^12-1.25;
    ch2 = 2.5*ch2/2^12-1.25;
else
    ch1 = readtable(filenameA);
    ch2 = readtable(filenameB);
    ch1 = table2array(ch1);
    ch2 = table2array(ch2);
end

ch1 = ch1 - mean(ch1);
ch2 = ch2 - mean(ch2);

Ch1 = fft(ch1);
Ch2 = fft(ch2);
%% Isolate and filter
L = length(ch1);
df = S/L; % frequency resolution
f = df*(0:(L-1)).'; % frequency points

% plot raw data and select windows
figure('units','normalized','outerposition',[0 0.05 1 0.95])
set(gcf,'color','w','DefaultAxesFontSize',24)
ax(1) = subplot(2,1,1);
semilogy(f(100:10:round(L/2))*10^-6,abs(Ch1(100:10:round(L/2))),'k')
axis tight
set(gca,'XTickLabel',[],'YTickLabel',[]);
ax(2) = subplot(2,1,2);
semilogy(f(100:10:round(L/2))*10^-6,abs(Ch2(100:10:round(L/2))),'k')
axis tight
set(gca,'YTickLabel',[]);
title(ax(1),'Select left and right bounds of dual-comb')
[fs(1:2), ~] = ginput(2);
Chs(1) = find(ismember(ax, gca));
title(ax(1),'Select left and right bounds of lowest freq CW beat note')
[fs(3:4), ~] = ginput(2);
Chs(2) = find(ismember(ax, gca));
answer1 = questdlg('That beat note was with which CW?','Identify','1','2','Cancel','Cancel');
title(ax(1),'Select left and right bounds of 2nd-lowest freq CW beat note')
[fs(5:6), ~] = ginput(2);
Chs(3) = find(ismember(ax, gca));
answer2 = questdlg('That beat note was with which CW?','Identify','1','2','Cancel','Cancel');
title(ax(1),'Select left and right bounds of 3rd-lowest freq CW beat note')
[fs(7:8), ~] = ginput(2);
Chs(4) = find(ismember(ax, gca));
answer3 = questdlg('That beat note was with which CW?','Identify','1','2','Cancel','Cancel');
title(ax(1),'Select left and right bounds of highest freq CW beat note')
[fs(9:10), ~] = ginput(2);
Chs(5) = find(ismember(ax, gca));
answer4 = questdlg('That beat note was with which CW?','Identify','1','2','Cancel','Cancel');
close
fs = 2*10^6*fs/S; % normalized freqs for iir

% Analytic BP filter
data_A = zeros(L,5);
for i = 1:5                              % 1 = DCS, 2 = low CW, 3 = high CW
    if fs(2*i) > 0.90
        [b, a] = butter(8,fs(2*i-1),'high');
    elseif fs(2*i-1) < 0.10
        [b, a] = butter(8,fs(2*i),'low');
    else
        [b, a] = butter(8,[fs(2*i-1) fs(2*i)],'bandpass');
    end
    if Chs(i) == 1
        data_A(:,i) = filtfilt(b,a,ch1);
    else
        data_A(:,i) = filtfilt(b,a,ch2);
    end
    data_A(:,i) = hilbert(data_A(:,i));
end

% Account for optical and electrical delays
data_A(:,1) = circshift(data_A(:,1),round(t_delay*S/2),1);
data_A(:,2) = circshift(data_A(:,2),round(-t_delay*S/2),1);
data_A(:,3) = circshift(data_A(:,3),round(-t_delay*S/2),1);
data_A(:,4) = circshift(data_A(:,4),round(-t_delay*S/2),1);
data_A(:,5) = circshift(data_A(:,5),round(-t_delay*S/2),1);
data_A = data_A(1+round(abs(t_delay)*S/2):end-round(abs(t_delay)*S/2),:);
if plots(1) == 1
    ch1 = ch1(1+round(abs(t_delay)*S/2):end-round(abs(t_delay)*S/2),:);
    ch2 = ch2(1+round(abs(t_delay)*S/2):end-round(abs(t_delay)*S/2),:);
    Ch1 = fft(ch1);
    Ch2 = fft(ch2);
end

% For nice decimation
while mod(length(data_A(:,1)),r) ~= 0
    data_A = data_A(1:end-1,:);
    ch1 = ch1(1:end-1);
    ch2 = ch2(1:end-1);
end

% update axes
L = length(data_A(:,1));
dt = 1/S; % time resolution
t = dt*(0:(L-1)).'; % time points
df = S/(L); % frequency resolution
f = df*(0:(L-1)).'; % frequency points
%% Plot filtered data
if plots(1) == 1
    figure('units','normalized','outerposition',[0 0.05 1 0.95])
    set(gcf,'color','w','DefaultAxesFontSize',24)
    
    Ch1 = fft(ch1);
    Ch2 = fft(ch2);
    
    subplot(2,2,1)
    plot(t*10^3,real(data_A(:,1)),'color',0.5*[1 0 1]) %,'color',[0.85 0.325 0.098])
    axis tight
    title('Filtered dual-comb')
    xlabel('t (ms)')
    ylabel('y (V)')
    
    subplot(2,2,2)
    hold on
    if Chs(1) == 1
        plot(f*10^-6,2*abs(Ch1),'k')
    else
        plot(f*10^-6,2*abs(Ch2),'k')
    end
    plot(f*10^-6,abs(fft(data_A(:,1))),'color',0.5*[1 0 1])
    hold off
    axis tight
    set(gca,'Ytick',[],'YColor','none')
    title('Filter Check - DCS')
    xlabel('f (MHz)')
    ylabel('y (a.u.)')

    subplot(2,1,2)
    hold on
    if Chs(2) == 1
        plot(f*10^-6,2*abs(Ch1),'k')
    else
        plot(f*10^-6,2*abs(Ch2),'k')
    end
    trans = 0.5;
    h = plot(f*10^-6,abs(fft(data_A(:,2))));
    col = h.Color;
    h.Color = [col, trans];
    h = plot(f*10^-6,abs(fft(data_A(:,3))));
    col = h.Color;
    h.Color = [col, trans];
    h = plot(f*10^-6,abs(fft(data_A(:,4))));
    col = h.Color;
    h.Color = [col, trans];
    h = plot(f*10^-6,abs(fft(data_A(:,5))));
    col = h.Color;
    h.Color = [col, trans];
    hold off
    axis tight
    set(gca,'Ytick',[],'YColor','none')
    title('Filter Check - CW beat notes')
    xlabel('f (MHz)')
    ylabel('y (a.u.)')
end
%% Mate the animals
% Remove amplitude noise from CW
data_A(:,2) = data_A(:,2)./abs(data_A(:,2));
data_A(:,3) = data_A(:,3)./abs(data_A(:,3));
data_A(:,4) = data_A(:,4)./abs(data_A(:,4));
data_A(:,5) = data_A(:,5)./abs(data_A(:,5));

% Assign the beat notes
cw1 = [];
cw2 = [];
if strcmp(answer1,'1')
    cw1 = [cw1, 1];
elseif strcmp(answer1,'2')
    cw2 = [cw2, 1];
end
if strcmp(answer2,'1')
    cw1 = [cw1, 2];
elseif strcmp(answer2,'2')
    cw2 = [cw2, 2];
end
if strcmp(answer3,'1')
    cw1 = [cw1, 3];
elseif strcmp(answer3,'2')
    cw2 = [cw2, 3];
end
if strcmp(answer4,'1')
    cw1 = [cw1, 4];
elseif strcmp(answer4,'2')
    cw2 = [cw2, 4];
end

% Assign red/blue
answer = questdlg('Which is the redder laser?','Colors','CW1','CW2','Cancel','CW1');
if strcmp(answer,'CW1')
    cwR = cw1 + 1;
    cwB = cw2 + 1;
elseif strcmp(answer,'CW2')
    cwR = cw2 + 1;
    cwB = cw1 + 1;
else
    error('Red/Blue bad')
end

% Mate beat notes to produce offspring
answerR = questdlg('Which directions do CW1 beat notes move?','Beat notes',...
    'Same','Opposite','Cancel','Same');
answerB = questdlg('Which directions do CW2 beat notes move?','Beat notes',...
    'Same','Opposite','Cancel','Same');

if strcmp(answerR,'Same')
    toothR_A = conj(data_A(:,cwR(1))).*data_A(:,cwR(2));
else
    toothR_A = data_A(:,cwR(1)).*data_A(:,cwR(2));
end
if strcmp(answerB,'Same')
    toothB_A = conj(data_A(:,cwB(1))).*data_A(:,cwB(2));
else
    toothB_A = data_A(:,cwB(1)).*data_A(:,cwB(2));
end

clearvars fs
% plot offspring and select window
figure('units','normalized','outerposition',[0 0.05 1 0.95])
set(gcf,'color','w','DefaultAxesFontSize',24)

subplot(2,1,1);
plot(f*10^-6,abs(fft(toothR_A)),'r')
axis tight
set(gca,'YTickLabel',[]);
title('Select left and right bounds of red offspring')
[fs(1:2), ~] = ginput(2);

subplot(2,1,2);
plot(f*10^-6,abs(fft(toothB_A)),'b')
axis tight
set(gca,'YTickLabel',[]);

title('Select left and right bounds of blue offspring')
[fs(3:4), ~] = ginput(2);
close
fs = 2*10^6*fs/S;

if plots(2) == 1
    toothR_A_old = toothR_A;
    toothB_A_old = toothB_A;
end

if fs(1) < 0.10
    [b, a] = butter(8,fs(2),'low');
else
    [b, a] = butter(8,[fs(1) fs(2)],'bandpass');
end
toothR_A = filtfilt(b,a,toothR_A);
if fs(3) < 0.10
    [b, a] = butter(8,fs(4),'low');
else
    [b, a] = butter(8,[fs(3) fs(4)],'bandpass');
end
toothB_A = filtfilt(b,a,toothB_A);

window = ones(L,1);
window(round(L/2)+1:end) = 0;
toothR_A = ifft(window.*fft(toothR_A));
toothB_A = ifft(window.*fft(toothB_A));
if plots(2) == 1
    toothR_A_new = toothR_A;
    toothB_A_new = toothB_A;
end
% Remove amplitude noise
toothR_A = toothR_A./abs(toothR_A);
toothB_A = toothB_A./abs(toothB_A);
%% Correct offset noise
L = length(data_A(:,1));
% refresh axes
dt = 1/S; % time resolution
t = dt*(0:(L-1)).'; % time points
df = S/L; % frequency resolution
f = df*(0:(L-1)).'; % frequency points

clock_CW = toothB_A.*conj(toothR_A);
dcs_A = data_A(:,1).*conj(toothR_A);
%% Plot offspring and offset corrected signal
if plots(2) == 1
    figure('units','normalized','outerposition',[0 0.05 1 0.95])
    set(gcf,'color','w','DefaultAxesFontSize',24)
    subplot(2,2,1)
    hold on
    plot(f*10^-6,abs(fft(toothR_A_old)),'k')
    plot(f*10^-6,abs(fft(toothR_A_new)),'r')
    hold off
    axis tight
    set(gca,'Ytick',[],'YColor','none')
    title('Filter Check - Red Offspring')
    xlabel('f (MHz)')
    ylabel('A (a.u.)')
    
    subplot(2,2,3)
    hold on
    plot(f*10^-6,abs(fft(toothB_A_old)),'k')
    plot(f*10^-6,abs(fft(toothB_A_new)),'b')
    hold off
    axis tight
    set(gca,'Ytick',[],'YColor','none')
    title('Filter Check - Blue Offspring')
    xlabel('f (MHz)')
    ylabel('A (a.u.)')
    
    subplot(2,2,2)
    plot((f-S/2)*10^-6,abs(fftshift(fft(data_A(:,1).*conj(toothR_A)))),'r')
    axis tight
    set(gca,'Ytick',[],'YColor','none')
    title('Corrected by red offspring')
    xlabel('f (MHz)')
    ylabel('A (a.u.)')
    
    subplot(2,2,4)
    plot((f-S/2)*10^-6,abs(fftshift(fft(data_A(:,1).*conj(toothB_A)))),'b')
    axis tight
    set(gca,'Ytick',[],'YColor','none')
    title('Corrected by blue offspring')
    xlabel('f (MHz)')
    ylabel('A (a.u.)')
end
%% Condition CW clock
% plot and select CW clock window
figure('units','normalized','outerposition',[0 0.05 1 0.95])
set(gcf,'color','w','DefaultAxesFontSize',24)
plot((f-S/2)*10^-6,abs(fftshift(fft(clock_CW))))       %Question: Why is f-S/2
axis tight
set(gca,'Ytick',[],'YColor','none')
% xlim([0 500]) % show first 10 harmonics
xlabel('f (MHz)')
title('Select left and right of CW clock')
[fs, ~] = ginput(2);
close
fs = 2*10^6*fs/S;

if plots(3) == 1
    clock_CW_old = clock_CW;
end

% isolate clock
[b, a] = butter(8,max(fs),'low');
clock_CW = filtfilt(b,a,clock_CW);
%% Plot clock
if plots(2) == 1
    figure('units','normalized','outerposition',[0 0.05 1 0.95])
    set(gcf,'color','w','DefaultAxesFontSize',24)
%     plot(f*10^-6,abs(fft(toothR_A)),'r')
    hold on
%     plot(f*10^-6,abs(fft(toothB_A)),'b')
    plot((f-S/2)*10^-6,abs(fftshift(fft(clock_CW_old))),'k')
    plot((f-S/2)*10^-6,abs(fftshift(fft(clock_CW))),'color',0.5*[1 0 1])
%     plot(f*10^-6,abs(fft(conj(clock2).*toothR_A)),'color',0.5*[0 1 1])
    hold off
    axis tight
    set(gca,'Ytick',[],'YColor','none')
    title('Filter check - Clock')
    xlabel('f (MHz)')
    ylabel('A (a.u.)')
end
%% Generate DDFG
ddfg = dcs_A.*conj(dcs_A);
ddfg = ddfg - mean(ddfg);
DDFG = fft(ddfg);

if plots(3) == 1
    DDFG_old = DDFG;
end
%% Unwrap clock
temp = fft(clock_CW);
temp = [temp(1:L/2); zeros(3*L,1); temp((L/2+1):end)];  % Question: what's the motivation of adding zero(3*L,1)
clock2_up = ifft(temp);
%% Determine ideal Delta_f
phase = unwrap(angle(clock2_up));
phase = decimate(phase,4);
m = 7;
Delta_f = (S/(2*pi))*((max(phase) - min(phase))/(L-1)); % exact average!!
% Now choose new Delta_f to force L_burst = integer IMPORTANT
L_burst = ceil(S/Delta_f);
Delta_f = S/L_burst;
%% Resample
N_bursts = floor(L/L_burst); % always drop the last one
if mod(N_bursts,2) ~= 0 % make sure it's even    Question: why do we need to make sure the number of bursts is even?
N_bursts = N_bursts - 1;
end
L = N_bursts*L_burst; % truncate record length IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% update axes
dt = 1/S; % time resolution
t = dt*(0:(L-1)).'; % time points
df = S/L; % frequency resolution
f = df*(0:(L-1)).'; % frequency points

phase_new = min(phase) + (2*pi/S)*Delta_f.*(0:(L-1)).'; % target phase ramp
inds = interp1(phase,(1:length(dcs_A)).',phase_new,'pchip'); % corresponding indices
dcs_A = interp1((1:length(dcs_A)).',dcs_A,inds,'pchip'); % resample & truncation

if plots(3) == 1
    % regenerate ddfg
    ddfg = dcs_A.*conj(dcs_A);
    ddfg = ddfg - mean(ddfg);
    DDFG = fft(ddfg);           % Question: Why do we need to regenerate DDFG here?
    clock_CW_new = interp1((1:length(clock_CW)).',clock_CW,inds,'pchip');
end
%% Plot resample verification
if plots(3) == 1
    figure('units','normalized','outerposition',[0 0.05 1 0.95])
    set(gcf,'color','w','DefaultAxesFontSize',24)
    subplot(1,2,1)
    plot(((S/length(clock_CW_old))*(0:(length(clock_CW_old)-1)).')*10^-6,2*abs(clock_CW_old),'k')
    hold on
    plot(f*10^-6,2*abs(clock_CW_new),'color',0.5*[1 0 1])
    hold off
    axis tight
    set(gca,'Ytick',[],'YColor','none')
    title('CW clock before and after')
    xlabel('f (MHz)')
    
    subplot(1,2,2)
    plot(((S/length(DDFG_old))*(0:(length(DDFG_old)-1)).')*10^-6,2*abs(DDFG_old),'k')
    hold on
    plot(f*10^-6,2*abs(DDFG),'color',0.5*[1 0 1])
    hold off
    axis tight
    set(gca,'Ytick',[],'YColor','none')
    title('DDFG before and after')
    xlabel('f (MHz)')
end
%% Pull peaks
[~, ind] = min(abs(f-Delta_f));
ind0 = L/2 + 1;
inds = ((ind-1):(ind-1):L/2).';
inds = [-flipud(inds); 0; inds];
inds = inds + ind0;
inds = inds(2:end-1);
Comb = abs(fftshift(fft(dcs_A)));
peaks = Comb(inds);
%% Plot fully corrected data
if plots(4) == 1
    figure('units','normalized','outerposition',[0 0.05 1 0.95])
    set(gcf,'color','w','DefaultAxesFontSize',24)
    subplot(1,2,1)
    plot(t*10^3,real(dcs_A),'color',0.5*[1 0 1])
    axis tight
    title('Dual-comb')
    xlabel('t (ms)')
    ylabel('y (V)')
    
    subplot(1,2,2)
    plot((f-S/2)*10^-6,Comb,'color',0.5*[1 0 1])
    hold on
    plot((f(inds)-S/2)*10^-6,peaks,'k')
    hold off
    axis tight
    yl = ylim;
    ylim([yl(2)/10^3 yl(2)])
    xlim(10^-6*Delta_f*([-30 30]))
    set(gca,'Ytick',[],'YColor','none')
    xlabel('f (MHz)')
    title(strcat('N_{bursts} = ',num2str(N_bursts),', \Delta_{f rep} = ',num2str(10^-6*round(Delta_f)),'... MHz'))
end
%% Save data
answer = questdlg('Save fully corrected data?','Save',...
    'Yes','No','No');

if strcmp(answer,'Yes')
    csvwrite(strcat(filenameA(1:end-4),'_processed.csv'),dcs_A)
    csvwrite(strcat(filenameA(1:end-4),'_peaks_only.csv'),peaks)
    export_fig(strcat(filenameA(1:end-4),'_plot.png'), '-r100')
end
