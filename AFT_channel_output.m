function r = AFT_channel_output(N_AFT, Num_OFDM_sym, taps, delay_taps, Doppler_taps, chan_coef,sigma_2,x_with_CP)
%% wireless channel and noise
N_CP = max(delay_taps);
s_chan_mat = zeros(size(x_with_CP));
for i = 1:Num_OFDM_sym
    for itao = 1:taps
        %         s_chan_mat(i, :) = s_chan_mat(i, :)+exp(1i*2*pi*Doppler_taps(itao)*delay_taps(itao))*chan_coef(itao)*exp(-1j*2*pi ...
        %             *(-N_CP:N_AFT-1)*Doppler_taps(itao)).*circshift(x_with_CP(i, :) ,delay_taps(itao), 2);
        s_chan_mat(i, :) = s_chan_mat(i, :)+chan_coef(itao)*exp(-1j*2*pi ...
            *(-N_CP:N_AFT-1)*Doppler_taps(itao)).*circshift(x_with_CP(i, :) ,delay_taps(itao), 2);
    end
end
noise = sqrt(sigma_2/2)*(randn(size(s_chan_mat)) + 1i*randn(size(s_chan_mat)));
r = s_chan_mat + noise;
r = r(:, N_CP+1:end);%discard cp
end