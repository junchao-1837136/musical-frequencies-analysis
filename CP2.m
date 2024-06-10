load 'CP2_SoundClip.mat'

Fs = 44100; %sample rate of the sound clip

S = y'; % transposes in order to have consistent dimensions
w = length(y)/4; % break the spectogram up into four time windows, otherwise it
% will be too big for MATLAB/autograder to run.

S1 = S((1-1)*w+1:1*w); % this will isolate the correct window
S2 = S((2-1)*w+1:2*w); % this will isolate the correct window
S3 = S((3-1)*w+1:3*w); % this will isolate the correct window
S4 = S((4-1)*w+1:4*w); % this will isolate the correct window

L = length(S1)/Fs; % length of each window in seconds
n = length(S1);  % number of elements in each window
t = [0:1/Fs:L - 1/Fs]; % t in sec. relative to the start of the window
tau = 0:0.1:L; % discretization for the Gabor transform
k = 2*pi*(1/L/2)*[0:n/2-1 -n/2:-1]; % discretization in frequency space
ks = fftshift(k); % gotta shift them freqs.

Sgt_spec = zeros(length(ks),length(tau)); % initializing the function for the spectrogram

%Gabor Transform Parameters
a = 400; %this will give you the correct width, so exp(-a(...))
range = [1:1800]; %use this when finding the max and index of the transformed Gabor filtered signal
% i.e., max(TransformedSignal_GaborFiltered(range))


% Repeat this part for S1, S2, S3, and S4.
% For this part we are going to make a spectrogram, but only for the freqs.
% of interest.  So first we'll do a Gabor transform, then we'll filter
% around our peak freq. in the regime of interest, and then we'll look at
% the spectrogram of that function.

% Gabor transform each S, just like we did in the lecture
% Week4_Spectrograms.m lines 141 to 146.
% You'll have to add code between line 144 and Sgt_spec at line 145.

% After line 144 find the index of your peak frequency (in absolute value)
% within the range of interest (this range is very forgiving so you don't
% have to match the autograder exactly, just use your judgement based on
% the figure in the assignment statement).  Then make a filter centered
% about the peak frequency.  Filter the Gabor transformed function.
% this is the function you will use in line 145 from the lecture to find
% your Sgt_spec.

for j = 1:length(tau)
   g = exp(-a*(t - tau(j)).^2); 
   Sg = g.*S1;
   Sgt = fft(Sg);
   [value, index] = max(abs(Sgt(range)));
   max_k = abs(k(index));
   filter = exp(-1/L*(abs(k)-max_k).^2);
   Sgt_spec(:,j) = fftshift(abs(filter.*Sgt)); 
end
ans1 = Sgt_spec;
for j = 1:length(tau)
   g = exp(-a*(t - tau(j)).^2); % Window function
   Sg = g.*S2;
   Sgt = fft(Sg);
   [value, index] = max(abs(Sgt(range)));
   max_k = abs(k(index));
   filter = exp(-1/L*(abs(k)-max_k).^2);
   Sgt_spec(:,j) = fftshift(abs(filter.*Sgt)); % We don't want to scale it
end
ans2 = Sgt_spec;
for j = 1:length(tau)
   g = exp(-a*(t - tau(j)).^2); % Window function
   Sg = g.*S3;
   Sgt = fft(Sg);
   [value, index] = max(abs(Sgt(range)));
   max_k = abs(k(index));
   filter = exp(-1/L*(abs(k)-max_k).^2);
   Sgt_spec(:,j) = fftshift(abs(filter.*Sgt)); % We don't want to scale it
end
ans3 = Sgt_spec;
for j = 1:length(tau)
   g = exp(-a*(t - tau(j)).^2); % Window function
   Sg = g.*S4;
   Sgt = fft(Sg);
   [value, index] = max(abs(Sgt(range)));
   max_k = abs(k(index));
   filter = exp(-1/L*(abs(k)-max_k).^2);
   Sgt_spec(:,j) = fftshift(abs(filter.*Sgt)); % We don't want to scale it
end
ans4 = Sgt_spec;



% Save Sgt_spec as variable A1 after your for loop.  Repeat for S2, S3, and
% S4 and save those Sgt_spec as A2, A3, and A4.  You don't have to rewrite
% the code, just copy and paste and use the respective S's, or write a for
% loop that iterates through S1 to S4.

A1 = ans1;   % Shape:  484560x110 double
A2 = ans2;   % Shape:  484560x110 double
A3 = ans3;   % Shape:  484560x110 double
A4 = ans4;   % Shape:  484560x110 double


% Plot of spectrogram for each window (for the report, not for autograder)
% just like we did in the lecture, but change your ylim to be in the range
% of interest for the sound clip.
figure();
pcolor(tau,ks,Sgt_plot)
shading interp
set(gca,'ylim',[0 700],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (sec.)'), ylabel('frequency (k)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Isolate the bassline
S = y'; % transposes in order to have consistent dimensions
L = length(S)/Fs; % total length in sec. 
n = length(S); % total number of elements in S
t = [0:1/Fs:L - 1/Fs]; % time discretization
k = (1/L)*[0:n/2-1 -n/2:-1]; % freq. discretization

% Take the Fourier transform of S, and in freq. space isolate all freqs. 
% (in absolute value) that you determine should be part of the baseline 
% according to spectrogram (or also just by listening); that is, all points
% in the transformed function not within the frequency range you determined
% should be set to zero (kind of like a Shannon filter, but simpler than
% what we did in lecture).
% You may have to do this part a few times with different thresholds to get
% it right.

Int = fft(S);
filter = zeros(size(Int));
for j = 1:n
    if abs(k(j))>250
        filter(j)=0;
    else
        filter(j)=1;
    end
end

Intf = Int.*filter;
Inf = abs(ifft(Intf));


% After thresholding the transformed function, take the inverse transform
% of the thresholded function and save it as A5.

A5 = Inf';     %Shape:  1938240x1 double


%Play sound (not for autograder)


%Plot the amplitude S over time (for the report, not for the autograder)
figure();
plot(t, A5);
xlabel('Time (sec.)');
ylabel('Amplitude');
title('Amplitude of Bassline over Time');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Isolate the guitar

%Same exact process as the baseline above, but you'll have to be more
%careful about the frequency range.
S = y'; %reinitialize the S from the previous part above.

Int = fft(S);
filter = zeros(size(Int));
for j = 1:n
    if abs(k(j))>1200 || abs(k(j))<83
        filter(j)=0;
    else
        filter(j)=1;
    end
end

Intf = Int.*filter;
Inf = abs(ifft(Intf));


A6 = Inf';     %Shape:  1938240x1 double

%Play sound (not for autograder)


%Plot the amplitude S over time (for the report, not for the autograder)

figure();
plot(t, A6);
xlabel('Time (sec.)');
ylabel('Amplitude');
title('Amplitude of Guitar over Time');