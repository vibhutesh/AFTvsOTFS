function x_est = AFT_mp_detector(N_AFT, Num_OFDM_sym, c0, c1, c2,taps,delay_taps,Doppler_taps,chan_coef, y)
N_CP        = max(delay_taps);
V_n_l       = zeros(N_AFT + N_CP);
G_n_l       = zeros(N_AFT + N_CP);
G_mat       = zeros(N_AFT);
x_est       = zeros(Num_OFDM_sym, N_AFT);
%H_Channel_eq   = zeros(N_AFT);
for n = -N_CP:N_AFT-1
    for l = 0:taps-1%Number of Cluster
        A_i_n                           = chan_coef(l+1);% It can be RV
        f_i_n                           = Doppler_taps(l+1);
        G_n_l(n+N_CP+1, delay_taps(l+1)+1)   = A_i_n*exp(-1i*2*pi*f_i_n*n);
        V_n_l(n+N_CP+1, delay_taps(l+1)+1) = G_n_l(n+N_CP+1, delay_taps(l+1)+1)*exp(1i*2*pi*c1*(delay_taps(l+1)^2 - 2*n*delay_taps(l+1)) + 1i*2*pi*(c0*n));
        if n >=0
            G_mat(n+1 , mod(n-delay_taps(l+1), N_AFT) + 1) = V_n_l(n+N_CP+1, delay_taps(l+1)+1);
            %Coeff  = 1;
            %if ((n <= (N_CP)) && (mod(n-delay(l+1), N) >= N - N_CP))
            %   Coeff  = exp(-1i*2*pi*c1*(N^2 - 2*N*(N - mod(n-delay(l+1), N))));
            %end
            %H_Channel_eq(n+1 , mod(n-delay(l+1), N) + 1) = G_n_l(n+N_CP+1, delay(l+1)+1)*Coeff;
        end
    end
end
F                   = 1/sqrt(N_AFT)*dftmtx(N_AFT);
Gamma_c2 = eye(N_AFT); % because c2 = 0;
H_eq                                                  = Gamma_c2*F*G_mat*F'*Gamma_c2';
h_channel_eq_total                                    = transpose(diag(H_eq));
h_channel_eq_total(ismember(h_channel_eq_total , 0))  = .0001;
for i=1:Num_OFDM_sym
    x_est(i, :) = y(i, :)./(h_channel_eq_total);
end
end