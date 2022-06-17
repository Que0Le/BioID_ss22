data = load("DavisData.mat");

%% 1.a
mean_uf = mean(data.xF);
mean_um = mean(data.xM);

%% 1.b
cov_sig_F = cov(data.xF)
cov_sig_M = cov(data.xM);

%% 1.c
var_Fm = cov_sig_F(1);
var_Fh = cov_sig_F(4);
var_Mm = cov_sig_M(1);
var_Mh = cov_sig_M(4);
r_Fmh = (cov_sig_F(2))/sqrt((cov_sig_F(1))*(cov_sig_F(4)))
r_Mmh = (cov_sig_M(2))/sqrt((cov_sig_M(1))*(cov_sig_M(4)));
% p_F = corrcoef(data.xF)
% p_M = corrcoef(data.xM);

%% 1.d
% mean vectors: hold the average values from the data sets (sum of all values/number of
% values). With these mean values of student's weight and hight we know
% that a student (Female or male) has an average hight/weight of x cm/kg. This
% value is not neccessary a real value from the data set!
%
% cov: according to wiki, covariance is a measure of the joint variability
% of two random variables. In our examle, these 2 variables are weight and
% height and the calculated variance values are located in the covariance 
% matrix 2x2. The (1,2) and (2,1) variance values are the same and tell us how 
% related to each other are the height and weight value of a student, while
% (1,1) and (2,2) means how the height or weight values are distributed
% compared to the mean value. In this data set, the weight and height are
% related to each other and self in the way that the higher a person is,
% the heavier he/she gets (represented by the positive variance values).  
% 
% correlation coefficients: These values (F: 0.513729293275502, M: 0.539290613882818)
% tell us how "good" are the progressive relations between height and
% weight. Female student's weight tends to be less increase with height.
%
%



%% 2.

x1 = [62, 168];
x2 = [64, 170];
x3 = [66, 172];
x4 = [68, 174];
p1_x1_F = mvnpdf(x1, mean_uf, cov_sig_F)
p1_x1_M = mvnpdf(x1, mean_um, cov_sig_M)

p2_x2_F = mvnpdf(x2, mean_uf, cov_sig_F)
p2_x2_M = mvnpdf(x2, mean_um, cov_sig_M)

p3_x3_F = mvnpdf(x3, mean_uf, cov_sig_F)
p3_x3_M = mvnpdf(x3, mean_um, cov_sig_M)

p4_x4_F = mvnpdf(x4, mean_uf, cov_sig_F)
p4_x4_M = mvnpdf(x4, mean_um, cov_sig_M)

%% 2.a equal priors 𝑃 𝜔1 = 𝑃 𝜔2:
% p1_x1_F = 0.0035 > p1_x1_M = 6.8266e-04 => x1 F
% p2_x2_F = 0.0025 > p2_x2_M = 0.0011     => x2 F
% p3_x3_F = 0.0015 = p3_x3_M = 0.0015     => x3 M (random, since P=P and p=p)
% p4_x4_F = 8.2584e-04 < p4_x4_M = 0.0019 => x4 M


%% 2.b prior probabilities of 𝑃(𝜔𝐹) = 0.4 and 𝑃(𝜔𝑀) = 0.6
P_F = 0.4;
P_M = 0.6;
pP1 = p1_x1_F * P_F - p1_x1_M * P_M;
pP2 = p2_x2_F * P_F - p2_x2_M * P_M;
pP3 = p3_x3_F * P_F - p3_x3_M * P_M;
pP4 = p4_x4_F * P_F - p4_x4_M * P_M;
% pP1 = 0.0010      > 0 => x1 F
% pP2 = 3.6882e-04  > 0 => x2 F
% pP3 = -2.7016e-04 < 0 => x3 M
% pP4 = -8.0005e-04 < 0 => x4 M


%% 2.c Assuming (reasonably) that 𝜆21 > 𝜆11 and 𝜆11 = 𝜆22 = 0
lambda_mf = 0.75;
lambda_fm = 1.5;
% p(𝒙|𝜔1)/p(𝒙|𝜔2) > (𝜆12 𝑃 𝜔2)/(𝜆21 𝑃 𝜔1) then decide 𝜔1 else decide 𝜔2
% x1: two_c_class_1 = 2.1661  > 0        => F
two_c_class_1 = p1_x1_F/p1_x1_M - (lambda_fm*P_M)/(lambda_mf*P_F)
% x2: two_c_class_2 = -0.6256 < 0        => M
two_c_class_2 = p2_x2_F/p2_x2_M - (lambda_fm*P_M)/(lambda_mf*P_F)
% x3: two_c_class_3 = -1.9566 < 0        => M
two_c_class_3 = p3_x3_F/p3_x3_M - (lambda_fm*P_M)/(lambda_mf*P_F)
% x4: two_c_class_4 = -2.5617 < 0        => M
two_c_class_4 = p4_x4_F/p4_x4_M - (lambda_fm*P_M)/(lambda_mf*P_F)


