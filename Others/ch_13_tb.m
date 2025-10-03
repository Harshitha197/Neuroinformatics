srate = 500;
f1     = [2,30];
time  = -1:1/srate:1; 
wavelet = exp(2*pi*1i*f.*time) .* exp(-time.^2./(2*s^2)); 

figure
subplot(221)
plot3(time,real(wavelet),imag(wavelet),'m')
xlabel('Time (ms)'), ylabel('real axis')
view(0,90)
title('Projection onto real and time axes')

subplot(222)
plot3(time,real(wavelet),imag(wavelet),'g')
xlabel('Time (ms)'), ylabel('imaginary axis')
view(0,0)
title('Projection onto imaginary and time axes') 
 
subplot(223)
plot3(time,real(wavelet),imag(wavelet),'k')
ylabel('real axis'), zlabel('imag axis')
view(90,0)
title('Projection onto imaginary and time axes')

subplot(224)
plot(time,real(wavelet),'b')
hold on
plot(time,imag(wavelet),'b:')
legend({'real part';'imaginary part'})
