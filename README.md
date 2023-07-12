# Ambiguity-Function
# implementation of ambiguity function in Matlab language and analyses of some signals
Ambiguity function is a 2D function of delay and Doppler frequency. It is an integral like fft integral. For narrowband
siganls we can use FFT algorithm instead of computing the integral; It is very faster because of its algorithm order(NlogN vs N^2).
When there is radial velocity between radar and target, we have frequency shift(Doppler frequency) in reflected signal.
In order to describe the response of the matched filter in the presence of Doppler frequency shift, the ambiguity function is used.
In fact, ambiguity funciton is the output of matched filter.
PS: To maximize detection possibility, we must maximize SNR. To maximize SNR we can use matched filter; So an important section
in designing all radar recievers is matched filter.
