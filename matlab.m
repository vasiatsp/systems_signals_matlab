%ASKHSH ENA - ANALYSH SYMATOS FONHS

% %audio reording rec = audiorecorder();
% disp('recording'); 
% recordblocking(rec,3);
% disp('end recording');
% 
% play(rec); 
% v = getaudiodata(rec);
% audiowrite('vasiliki.wav', v, 8000);

%ASKHSH 1
fsl = 8000; 
vas = audioread('vasiliki.wav'); 
t = linspace(0, 3,3*fsl ); 
figure(1)
plot(t, vas);

%normalization
nom = normalize(vas, 'range', [-1,1]); 
figure(2);
plot (t, nom); 

%hamming window 
n = linspace(0,200,2000);
M = 200; 
w = 0.54-0.46*cos(2*pi*n/M); 

%energy 
h= conv(nom.*nom, w);
nom = normalize(nom, 'range', [min(h), max(h)]);
figure(3);
plot (t, nom);
hold on
plot(linspace(0, 3, length(h)), h);
hold off

%zero crossing rate
k =rectwin(200); 
zcr = conv(abs(diff(sign(vas))),k);
figure(4)
plot(zcr); 


%d
figure(5)
subplot(3,1,1); 
plot(vas); 
subplot(3,1,2);
plot(h); 
subplot(3,1,3);
plot(zcr); 

%window 50ms
starting_point = fsl*1.29;
ending_point = fsl*1.35;
vas_window = vas(starting_point:ending_point);
figure(6);
t = linspace(1.5 , 1.55 , length(vas_window));
plot(t , vas_window);
sound(vas_window);


%dft
dtf_vas_window = fft(vas_window, 1024);
f = linspace (0, 8000, length(dtf_vas_window));
figure(7);
plot(f, abs(dtf_vas_window));
figure(8);
plot(f, 20*log10(abs(dtf_vas_window)));

%---ASKHSH 2.1-SXEDIASH BATHIMERATOU

% %FOR N=2
% a2 = [3,0];
% b2 = [1,1,1]; 
% 
% %FOR N=5
% a5 = [6,0]; 
% b5 = [1,1,1,1,1,1];

%FOR N = 2
a2 = (1:3);
b2 = (1:3);
a2(:) = 0; 
a2(1) = 1; 
b2(:) = 1/3; 

figure(1)
freqz(b2,a2); 

[z2, p2, k2] = tf2zp (b2, a2);
figure(10);
zplane(z2, p2);

%FOR N=5
a5 = (1:6);
b5 = (1:6);
a5(:) = 0; 
a5(1) = 1; 
b5(:) = 1/6; 

figure(8);
freqz(b5, a5);

[z5, p5, k5] = tf2zp (b5, a5);
figure(11);
zplane(z5, p5);

%FOR N=20
a20 = (1:21); 
b20 = (1:21);  
a20(:) = 0;
a20(1) = 1; 
b20(:) = 1/21; 

figure(9);
freqz(b20, a20);

[z20, p20, k20] = tf2zp (b20, a20);
figure(12);
zplane(z20, p20);

%Butterworth filter 
[bw2, aw2] = butter(2, 0.5); 
[bw8, aw8] = butter(8, 0.5); 
%FOR N=2
figure(13);
freqz(bw2, aw2); 

%FOR N=8
figure(14);
freqz(bw8, aw8); 

%ASKHSH 2.2-SXEDIASH ZONOPERATOY 
poles =[0.51+0.68i; 0.51-0.68i];
zeros=[1.5;-0.7]; 
figure(15)
zplane(zeros, poles); 

K = 0.15; 
[b,a] = zp2tf(zeros, poles, 0.15); 
figure(16)
freqz(b,a); 
 
figure(17)
impz(b, a);

figure(18)
stepz(b,a);

new_poles1 = [0.57+0.76i,0.57-0.76i];
new_poles2 = [0.6+0.8i, 0.6-0.8i]; 
new_poles3 = [0.63+0.84i, 0.63-0.84i]; 

[new_b_1 ,new_a_1] = zp2tf(zeros ,new_poles1 , 0.15);
figure(19)
zplane(new_b_1, new_a_1); 
%freqz(new_poles1, new_a_1);
figure(20)
stepz(new_poles1, new_a_1); 

[new_b_2 ,new_a_2] = zp2tf(zeros ,new_poles2 , 0.15);
figure(21)
zplane(new_b_2, new_a_2); 
%freqz(new_poles2, new_a_2);
figure(22)
stepz(new_poles2, new_a_2); 

[new_b_3 ,new_a_3] = zp2tf(zeros , new_poles3 , 0.15);
figure(23)
zplane(new_b_3, new_a_3); 
%freqz(new_poles3, new_a_3);
figure(24)
stepz(new_poles3, new_a_3); 

e_poles = [-0.68 + 0.51i; -0.68 - 0.51i];
[b_new , a_new] = zp2tf(zeros , e_poles , 0.15);
figure(25)
zplane(zeros, e_poles); 
figure(26)
freqz(b_new, a_new); 

%ASKHSH 3.1-ANALYSH MOUSIKVN FILTRWN
fsl=1600; %sampling frequency
Tf= 1/fsl; %sampling period 

flute = audioread('flute_note.wav');
brass = audioread("brass_note.wav");
string = audioread('string_note.wav'); 

%sound (flute);
%sound (brass); 
%sound (string); 

%original sygnal 
figure(251)
plot(linspace(0, fsl, length(flute)), flute); 
figure(261)
plot(linspace(0, fsl, length(brass)), brass); 
figure(271)
plot(linspace(0, fsl, length(string)), string);
 
%fft 
f_flute = fft(flute);
f_brass = fft(brass); 
f_string = fft(string);
figure(252)
plot(linspace(0, fsl, length(f_flute)), abs(f_flute));
figure(262)
plot(linspace(0, fsl, length(f_brass)), abs(f_brass));
figure(272)
plot(linspace(0, fsl, length(f_string)), abs(f_string));

brass_noisy = audioread('brass_noisy.wav'); 
f_brass_noisy = fft(brass_noisy);
%sound(abs(f_brass_noisy)); 
figure (1111)
plot(linspace(0, fsl, length(brass_noisy)), brass_noisy); 

%noise reduction (e)
f_cut= 0.3; % cutoff frequency
Wc = f_cut/(fsl/2); %frequency angle
[bw2, aw2] = butter(2, Wc); % exei pio aplh apokrish 
%[bw12, aw12] = butter(12, Wc);%to filtro taxhs 12 einai pio apotelesmatiko sthn apokoph ths syxnothtas 

brass_noisy_filter2 = filter(bw2, aw2, brass_noisy);
%brass_noisy_filter12 = filter(bw12, aw12, brass_noisy);

figure(28)
plot(linspace(0, fsl, length(brass_noisy_filter2)), brass_noisy_filter2);
sound(brass_noisy_filter2); 
% figure(29)
% plot(linspace(0, fsl, length(brass_noisy_filter12)), brass_noisy_filter12);
% %sound(abs(brass_noisy_filter12));

%bandpass filter
%frequency of 3th harmonic h = 590 Ηz previous next +- 190hz
h3 = 590; 
[n3,wn3] = buttord([h3-100 h3+100]/8000,[h3-500,h3+500]/8000,10,60);
[b_but3,a_but3]=butter(n3,wn3);
fillter3 = filter(b_but3,a_but3,string);
% frequency of th harmonic h = 980 Hz previous next +-190Hz
h5 = 980; 
[n5,wn5] = buttord([h5-100 h5+100]/8000,[h5-500 h5+500]/8000,10,60); 
[b_but5,a_but5]=butter(n5,wn5); 
filtter3 = filter(b_but3,a_but3,string);
fillter5 = filter(b_but5,a_but5,string);

figure(221)
subplot 211
plot (fillter3)
subplot 212
plot(fillter5)


%ASKHSH 3.2 SUNTHESH MOUSIKVN SHMATWN

flute = audioread('flute_note.wav'); 
ft = 200; 
len = 10*ft; 

% 
%apomonosh tmimatos eiso me 10 
% % seg_lenght = round((10*ft)*fsl); 
% seg2 = flute(flute :seg_lenght); 
% % figure(2)
% % plot(seg2); 
% seg_lenght = 10*ft;
% % seg2 = reshape(flute:seg_lenght), [], seg_lenght/ft); 
% % figure(8)
% plot(seg2); 

seg_lenght= 10*ft; 
seg = flute(1:seg_lenght); 
figure(29)
plot(seg)
n1 = 1000; 
n2 = 2000; 
samp = flute(n1:n2); 
figure
plot(samp)


fft_flute2 = fft(seg);
fft2 = abs(fft_flute2); 
freq = linspace(0, fsl, length(fft2)); 
figure(30)
plot(freq, fft2);
c_n =( 2*(fft2))/samp; 
figure (31)
stem(freq, c_n); 


% %phi calculation 
threshold = 100; 
l = length(fft_flute2); 
sig_harmonics = find(abs(fft_flute2) > threshold); 
phi = angle(fft_flute2(sig_harmonics))/l; 
phi_abjust = mod(phi + pi, 2*pi) - pi; 


%prosthesh hmitinoitdvn 
t = linspace(0, 0.16875, fsl); 
clear zeros ; 
x = 0; %anakataskeuasmeno shma 
    for n = 1:length(sig_harmonics)
        x = 0; 
        cn = c_n(n); %suntelesths platous
        p = phi(n); %phase
        fn = freq(n); %suxnothta
            x = x+cn*cos(2*pi*fn*t+p);
    end
    


figure(32)
plot( x); 
%sound(x,fsl)

audiowrite('flute.reconstructed.wav',x,fsl)
