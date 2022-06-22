data = load("SVdata.mat");
data.utterance;
titles = ["M1SET0", "M1SET1", "M3SET0", "M3SET1", "X1S1T0"];
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
fprintf("\n1.a ___________________________________________________________\n\n");
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

    fprintf('%s Log Short Time Energy: [%s]\n', titles(signal_th), join(string(stle_E(signal_th, 1:5)), ', '));
    fprintf('%s Voicing parameter    : [%s]\n', titles(signal_th), join(string(vp_VC(signal_th, 1:5)), ', '));
end


%% 1.b
fprintf("\n1.b ___________________________________________________________\n\n");

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
    fprintf('%s Number of vowel frames: %d\n', titles(signal_th), counter);
end


%% 1.c
fprintf("\n1.c ___________________________________________________________\n\n");

hW = hamming(frameWdth);
ceps_storage = {[], [], [], [], []};

for signal_th=1:5
    ceps = zeros(len_vowel_frames(signal_th, 1), 12);
    for v_i = 1:len_vowel_frames(signal_th, 1)

        start = (vowel_frames(signal_th, v_i)-1)*frameShft + 1;
        [y,ym] = rceps(hW.*data.utterance{1,signal_th}(start:(start-1+frameWdth)));
        %fprintf('a: [%s]\n', join(string(y(2:13)), ','));
        ceps(v_i, 1:12) = (y(2:13));
    end
    ceps_storage{1, signal_th} = ceps;
%     fprintf('%s 12D cepstral vector: [%s]\n', titles(signal_th), join(string(ceps(1:12, signal_th)), ', '));
end


%% 2.a
fprintf("\n2.a ___________________________________________________________\n\n");

% Concate ceps values from 2 records into 1.
ceps_M1 = vertcat(ceps_storage{1,1}, ceps_storage{1,2});
ceps_M3 = vertcat(ceps_storage{1,3}, ceps_storage{1,4});
ceps_X1 = ceps_storage{1,5};

mu_M1 = mean(ceps_M1);
mu_M3 = mean(ceps_M3);

sigma_M1 = transpose(diag(cov(ceps_M1)));
sigma_M3 = transpose(diag(cov(ceps_M3)));

fprintf('M1 mu: [%s]\n', join(string(mu_M1), ', '));
fprintf('M3 mu: [%s]\n', join(string(mu_M3), ', '));
fprintf('M1 sigma: [%s]\n', join(string(sigma_M1), ', '));
fprintf('M3 sigma: [%s]\n', join(string(sigma_M3), ', '));


%% 2.b
fprintf("\n2.b ___________________________________________________________\n\n");

log_likes_X1_M1 = zeros(1,10);
log_likes_X1_M3 = zeros(1,10);

for i=1:height(ceps_X1)

    log_likes_X1_M1(i) = log(mvnpdf(ceps_X1(i, 1:12), mu_M1, sigma_M1));
    log_likes_X1_M3(i) = log(mvnpdf(ceps_X1(i, 1:12), mu_M3, sigma_M3));

    fprintf('Log Likelihood for %dth frame in X: ', vowel_frames(5,i));
    fprintf('M1: %.4f dB  ', log_likes_X1_M1(i));
    fprintf('M2: %.4f dB\n', log_likes_X1_M3(i));
end


%% 2.c 
fprintf("\n2.c ___________________________________________________________\n\n");

mean_log_likes_X1_M1 = mean(log_likes_X1_M1);
mean_log_likes_X1_M3 = mean(log_likes_X1_M3);
LLR = mean_log_likes_X1_M1 - mean_log_likes_X1_M3;
LR = exp(LLR);
fprintf('Mean log Likelihood X1_M1: %.4f dB\n', mean_log_likes_X1_M1);
fprintf('Mean log Likelihood X1_M3: %.4f dB\n', mean_log_likes_X1_M3);
fprintf('LLR: %.4f\n', LLR);
fprintf('LR: %.4f\n', LR);





























