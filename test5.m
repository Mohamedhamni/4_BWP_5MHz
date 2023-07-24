fs = 7.68;
LTE_5_3_1.waveform = var.waveform;
tdw = LTE_5_3_1.waveform;
ps = PS5_768;

% Spectrum of the time domain waveform
figure
plot(linspace(-fs/2,fs/2,length(tdw)),20*log10(abs(fftshift(fft(tdw)))))
xlabel('Frequency (MHz)')
ylabel('Magnitude (dB)')

%% Balancing the array
ps_out_1 = conv(tdw,ps);
ps_out = ps_out_1(43:end-42);

figure
plot(linspace(-fs/2,fs/2,length(ps_out)),20*log10(abs(fftshift(fft(ps_out)))))
xlabel('Frequency (MHz)')
ylabel('Magnitude (dB)')

%% Upsampling
fs = 491.52;
ups_out = interp(ps_out,64);

figure
plot(linspace(-fs/2,fs/2,length(ups_out)),20*log10(abs(fftshift(fft(ups_out)))))
xlabel('Frequency (MHz)')
ylabel('Magnitude (dB)')


%% Shifting  for 4 Bandwidth part NCO operation

fcarr = 20;
fcarr2 = 40;
fcarr3 = 60;
fcarr4 = 80;
for ii = 1:length(ups_out)
   carr1_shifted(ii) = ups_out(ii).*exp(1i*2*pi*fcarr/fs*ii);
   carr2_shifted(ii) = ups_out(ii).*exp(1i*2*pi*fcarr2/fs*ii);
   carr3_shifted(ii) = ups_out(ii).*exp(1i*2*pi*fcarr3/fs*ii);
   carr4_shifted(ii) = ups_out(ii).*exp(1i*2*pi*fcarr4/fs*ii);
end
 carr_shifted = carr1_shifted + carr2_shifted + carr3_shifted + carr4_shifted;
 
 
figure
plot(linspace(-fs/2,fs/2,length( carr1_shifted)),20*log10(abs(fftshift(fft( carr1_shifted)))))
hold on;
title('Upsampled Waveform ');
%title('Shifted Waveform ');
plot(linspace(-fs/2,fs/2,length( carr2_shifted)),20*log10(abs(fftshift(fft( carr2_shifted)))))
hold on;
plot(linspace(-fs/2,fs/2,length( carr3_shifted)),20*log10(abs(fftshift(fft( carr3_shifted)))))
hold on;
plot(linspace(-fs/2,fs/2,length( carr4_shifted)),20*log10(abs(fftshift(fft( carr4_shifted)))))
xlabel('Frequency (MHz)')
ylabel('Magnitude (dB)')

figure
plot(linspace(-fs/2,fs/2,length( carr_shifted)),20*log10(abs(fftshift(fft( carr_shifted)))))
xlabel('Frequency (MHz)')
ylabel('Magnitude (dB)')

%% Cylic prifix removal of 4 BWP together


        td_sig = carr_shifted;
        Nsf = 10; % Number of subframes
        Nsym = 14*Nsf; % Number of symbols (14 symbols per subframe)
                Nsc = 300;
                CPint = 7;
                nFFT = 512*64; % FFT size
                cp0 = 40*64; % First CP length
                cpx = 36*64; % Normal CP length

        idx = 1;
        for sym_idx = 1:Nsym/7
            for sc_idx = 1:7
                if sc_idx == 1
                    cp_len = cp0;  % checking long                  
                else
                    cp_len = cpx;  % checking short
                end
                sym_grid(:,(sym_idx-1)*7+sc_idx) = td_sig(idx+cp_len:idx+cp_len+nFFT-1);
                idx = idx+cp_len+nFFT;
            end
        end
 cp_output = sym_grid;
 
 %% Appling fft for each resourse BLOCK
 
 % Input signal
 FFT_output = zeros(32768,140);
 for fftcol = 1:140
signal = cp_output(:,fftcol);
a = fft(signal);
% Apply FFT
FFT_output(:,fftcol) = a;
 end
 
 
%% Applig IFFT

% Calculation the range of each bandwidthpart
 
center_fft_1 = ceil(fcarr *1000/15);
center_fft_2 = ceil(fcarr2 *1000/15);
center_fft_3 = ceil(fcarr3 *1000/15);
center_fft_4 = ceil(fcarr4 *1000/15);

bandwidth_5 = ceil(5*1000/15);
% retreving the BWPs
% DC carrire check
BWP1 = FFT_output(center_fft_1-150: center_fft_1+ 255,:);
BWP2 = FFT_output(center_fft_2-256: center_fft_2+ 255,:);
BWP3 = FFT_output(center_fft_3-256: center_fft_3+ 255,:);
BWP4 = FFT_output(center_fft_4-256: center_fft_4+ 255,:);

% Appling IFFT

 IFFT_output_1 = zeros(512,140);
 for fftcol = 1:140
signal = BWP1(:,fftcol);
a = ifft(signal);
% Apply FFT
IFFT_output_1(:,fftcol) = a;
 end

  IFFT_output_2 = zeros(512,140);
 for fftcol = 1:140
signal = BWP2(:,fftcol);
a = ifft(signal);
% Apply FFT
IFFT_output_2(:,fftcol) = a;
 end
 
  IFFT_output_3 = zeros(512,140);
 for fftcol = 1:140
signal = BWP3(:,fftcol);
a = ifft(signal);
% Apply FFT
IFFT_output_3(:,fftcol) = a;
 end
 
  IFFT_output_4 = zeros(512,140);
 for fftcol = 1:140
signal = BWP4(:,fftcol);
a = ifft(signal);
% Apply FFT
IFFT_output_4(:,fftcol) = a;
 end








