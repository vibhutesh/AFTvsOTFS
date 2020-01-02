%in the name of Allah
clc
clear all
close all
tic
N                   = 4;%number of subchannel
NFFT                =N;
%
W                   = zeros(N);
a                   = 1;
b                   = 1;
d                   = 2;
c                   = (a*d - 1)/b;
k1                  = d/(4*pi*b);
k2                  = a/(4*pi*b);
Beta                = 1/(2*pi*b);
T                   = 1;
F                   = 1/(Beta*N*T);
delay               = [0 1];
power_channel_dB    = [0 0];
amp_k               = [1 1];
doppler_spread      =[0 1];
K                   =1/sqrt(N);% 1/sqrt(N);
c1                  = k1*T;
c2                  = k2*F;
f0                  = doppler_spread(1);
f1                  = doppler_spread(2);
l0                  = delay(1);
l1                  = delay(2);
c0                  = (f0*l1 - f1*l0)/(l1 - l0)
c1                  = (f0-f1)/(2*(l1 - l0));
c2                  = 0;
Gamma_c1            = zeros(N);
Gamma_c2            = zeros(N);
Gamma_c0            = zeros(N);
F                   = 1/sqrt(N)*dftmtx(N);
for n = 0:N-1
    Gamma_c1(n+1, n+1) = exp(-1i*2*pi*(c1*n^2));
    Gamma_c2(n+1, n+1) = exp(-1i*2*pi*(c2*n^2));
    Gamma_c0(n+1, n+1) = exp(1i*2*pi*(c0*n));
    for m = 0:N-1
        W(n+1, m+1) =K*exp(-1i*2*pi*(c2*n^2 + n*m/N + c1*m^2)); % W = GammaC2*F*GammaC1
    end
end
N_CP                = max(delay);
N_ofdm_sym          = NFFT + N_CP;
Tsym                = .001;
Ts                  = Tsym/(NFFT+N_CP);
num_bits            =4;%we need 16 signal for modulation
M                   = 2^num_bits;
SNR_dB              = 0:5:30;
Num_OFDM_sym        = 1;
ave_last            = 100;
power_signal        = 0;
G_mat               = zeros(N);
%transmitter
bit2de              = zeros(1, N*Num_OFDM_sym);
x_with_CP           = zeros(Num_OFDM_sym, N_ofdm_sym);
y_without_cp        = zeros(Num_OFDM_sym, N);
X_received_p        = zeros(Num_OFDM_sym, N);
X_received          = zeros(Num_OFDM_sym, N);
BER                 = zeros(1, length(SNR_dB));



for k = 0:length(SNR_dB)
    num_error   =0;
    num_bit     =0;
    for BER_ave =1:ave_last
        bits    = randi([0 1] ,1 ,N*num_bits*Num_OFDM_sym );%number of bit we need for Num_OFDM_sym ofdm symbol
        j       = 0;
        for i=1:num_bits:length(bits)-num_bits+1
            j           =j+1;
            bit2de(j)   = bi2de(bits(i :i+num_bits-1));
        end
        X           = qammod(bit2de , M);
        mean_En     = mean(abs(abs(X).^2));
        X           = X/mean_En;
        for i=1:Num_OFDM_sym
            x                   = transpose(sqrt(N)*transpose(diag(Gamma_c1)').*ifft(transpose(diag(Gamma_c2)'.*(X((i-1)*N+1:i*(N)))), NFFT));
            v                   = N_CP:-1:1;
            x_with_CP(i , :)    = [x(NFFT-N_CP+1:NFFT).*exp(-1i*2*pi*c1.*(N^2 - 2*N*(v))) x];%adding CP
        end
        %parallel to serial
        x_serial    = reshape(transpose(x_with_CP) ,[1,size(x_with_CP,1)*size(x_with_CP,2)]);
        
        s_AFT = AFT_modulation(N,Num_OFDM_sym, N_CP, c1, c2, X);
        
        %g_n_l
        V_n_l       = zeros(N_ofdm_sym);
        G_n_l       = zeros(N_ofdm_sym);
        H_Channel_eq   = zeros(N);
        for n = -N_CP:N-1
            for l = 0:length(delay)-1%Number of Cluster
                A_i_n                           = amp_k(l+1);% It can be RV
                f_i_n                           = doppler_spread(l+1);
                G_n_l(n+N_CP+1, delay(l+1)+1)   = A_i_n*exp(-1i*2*pi*f_i_n*n);
                V_n_l(n+N_CP+1, delay(l+1)+1) = G_n_l(n+N_CP+1, delay(l+1)+1)*exp(1i*2*pi*c1*(delay(l+1)^2 - 2*n*delay(l+1)) + 1i*2*pi*(c0*n));
                if n >=0
                    G_mat(n+1 , mod(n-delay(l+1), N) + 1) = V_n_l(n+N_CP+1, delay(l+1)+1);
                    Coeff  = 1;
                    if ((n <= (N_CP)) && (mod(n-delay(l+1), N) >= N - N_CP))
                        Coeff  = exp(-1i*2*pi*c1*(N^2 - 2*N*(N - mod(n-delay(l+1), N))));
                    end
                    H_Channel_eq(n+1 , mod(n-delay(l+1), N) + 1) = G_n_l(n+N_CP+1, delay(l+1)+1)*Coeff;
                end
            end
        end
        H_eq       = Gamma_c2*F*G_mat*F'*Gamma_c2';
        if max(max(abs(Gamma_c0*Gamma_c1*H_Channel_eq*Gamma_c1' - G_mat))) > 1e-9
            ali = 1;
        end
        H_mat_conv = zeros(N_ofdm_sym);
        for n = -N_CP:N-1
            for l = 0:length(delay)-1
                H_mat_conv(n+N_CP+1, mod(n+N_CP-delay(l+1), N_ofdm_sym) + 1) = G_n_l(n+N_CP+1, delay(l+1)+1);
            end
        end
        yy      = zeros(N_ofdm_sym, Num_OFDM_sym);
        for i = 1:Num_OFDM_sym
            yy(:, i) =  H_mat_conv*transpose(x_with_CP(i, :));
        end
        
        s_chan_mat = zeros(size(x_with_CP));
        for i = 1:Num_OFDM_sym
            for itao = 1:length(delay)
                s_chan_mat(i, :) = s_chan_mat(i, :)+amp_k(itao)*exp(-1j*2*pi ...
                    *(-N_CP:N-1)*doppler_spread(itao)).*circshift(x_with_CP(i, :) ,delay(itao));
            end
        end
        
        y       = reshape(yy, 1, []);
        if k==0
            power_signal = power_signal + sum(abs(y(1:Num_OFDM_sym*N_ofdm_sym)).^2);
            continue;
        end
        SNR         = 10^(SNR_dB(k)/10);
        noise       = 0;%sqrt(power_signal/(2*SNR))*(randn([1,length(y)])+1j*randn([1,length(y)]));
        y           = y+noise;
        h_channel_eq_total                           = transpose(diag(H_eq));
        h_channel_eq_total(ismember(h_channel_eq_total , 0))  = .0001;
        r_AFT = AFT_channel_output(N, Num_OFDM_sym, 2, delay, doppler_spread, amp_k,0,s_AFT);
        for i=1:Num_OFDM_sym
            y_without_cp(i , :)     = y((i-1)*N_ofdm_sym+N_CP+1:i*N_ofdm_sym);
            Y                       = (1/sqrt(N))*transpose((diag(Gamma_c2)).*fft(transpose(transpose(diag(Gamma_c1*Gamma_c0)).*(y_without_cp(i, :))), NFFT));
            X_received_p(i, :)      = Y./(h_channel_eq_total);
            X_received(i, :)        = X_received_p(i, 1:N);
        end
        x_serial_received           = reshape(transpose(X_received) ,[1,size(X_received,1)*size(X_received,2)]);
        bit_mat                     = de2bi(qamdemod(x_serial_received*mean_En , M), num_bits);
        bits_r                      = reshape(transpose(bit_mat) ,[1 , size(bit_mat,1)*size(bit_mat,2)]);
        num_error                   = num_error + sum(xor(bits_r , bits));
        num_bit                     = num_bit + length(bits_r);
    end
    if k == 0
        power_signal = power_signal/(Num_OFDM_sym*N_ofdm_sym*ave_last);
    else
        BER(k) = num_error/num_bit;
    end
end
BER
semilogy(SNR_dB ,BER )
hold on

toc


