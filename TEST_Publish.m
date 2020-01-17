%%
% This is a test

A1 = openfig('OTFSvsAFT_4QAM_N_8_M8_both.fig')
print(A1,'-dpdf')

A1 = openfig('OTFSvsAFT_4QAM_N_16_M16_MMSE.fig')
print(A1,'-dpdf')

A2 = openfig('OTFSvsAFT_16QAM_N_8_M_8_both.fig')
print(A2,'-dpdf')

A3 = openfig('OTFSvsAFT_16QAM_N_16_M16_MMSE.fig')
print(A3,'-dpdf')

A4 = openfig('OTFSvsAFT_64QAM_N_8_M8_MMSE.fig')
print(A4,'-dpdf')

A5 = openfig('OTFSvsAFT_64QAM_N_16_M16_MMSE.fig')
print(A5,'-dpdf')
