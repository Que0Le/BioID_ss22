data = load("SVdata.mat");
data.utterance;
% M1SET0 and M1SET1 spoken by a known speaker M1;
% M3SET0 and M3SET1 spoken by a known speaker M3
% X1S1T0 spoken by an unknown speaker X1

fsampl = 12500;
% Frame width [samples]: 
frameWdth=256; 
% Frame shift [samples]: 
frameShft=128; 
% Energy threshold [dB]: 
eThr=-45; 
% Voicing threshold: 
vcThr=0.4;
% Minimum fundamental frequency [Hz]: 
f0Min=80; 
% Maximum fundamental frequency [Hz]:
f0Max=200; 

%% 1.a
stle_E  = zeros(5,1);  % short-time log energy
len_stle_E = zeros(5,1);
vp_VC   = zeros(5,1);  % voicing parameter
len_vp_VC   = zeros(5,1);
kmin = floor(fsampl/f0Max);
kmax = ceil(fsampl/f0Min);
for signal_th = 1:5
    nbr_shifts = length(data.utterance{1,signal_th})/frameShft - 1;
    for i = 0:(nbr_shifts-1)
        sum = 0;
        xcorr_result = [];   % sd

        % Calculate short-time log energy 𝐸𝑖 in dB
        for cell_th = 1:frameWdth
            sum = sum + (data.utterance{1,signal_th}(i*frameShft + cell_th))^2;
        end
        stle_E(signal_th, i+1)= 10*log10(sum);

        % Calculate voicing param
        xcorr_result = xcorr(data.utterance{1,signal_th}((i*frameShft+1):(i*frameShft+frameWdth)),kmax,'normalized');
        vp_VC(signal_th, i+1)= max(xcorr_result(220:315));
    end
    len_stle_E(signal_th, 1) = nbr_shifts;
    len_vp_VC(signal_th, 1) = nbr_shifts;
end


%% 1.b
vowel_frames = zeros(5,1);   % hold the index of vowel frames
len_vowel_frames = zeros(5,1);
for signal_th=1:5
    counter = 0;
    for i = 1:len_stle_E(signal_th, 1)
        if stle_E(signal_th,i) > eThr && vp_VC(signal_th,i) > vcThr 
            counter = counter + 1;
            vowel_frames(signal_th, counter) = i;
        end
        len_vowel_frames(signal_th, 1) = counter;
    end
%     fprintf('Dataset %d \n', signal_th)
%     fprintf('Vowel Frames: %d\n',vowel_frames(signal_th));
end


%% 1.c
hW = hamming(frameWdth);
storage = {[], [], [], [], []};

for signal_th=1:5
    ceps = zeros(12, len_vowel_frames(signal_th, 1));
    for v_i = 1:len_vowel_frames(signal_th, 1)

        start = (vowel_frames(v_i)-1)*frameShft + 1;
        k = data.utterance{1,1}(start:(start-1+frameWdth));
        [y,ym] = rceps(hW.*data.utterance{1,signal_th}(start:(start-1+frameWdth)));
        ceps(1:12, v_i) = y(2:13);
    end
    storage{1, signal_th} = ceps;
end

