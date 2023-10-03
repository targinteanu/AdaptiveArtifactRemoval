function [wP, P, w, Y] = PowerSpectrum(y, Fs)
    % y: a row vector/matrix in time domain 
    % wP: one-sided frequency 
    % P: power spectrum (1-sided) of Y
    % w: two-sided frequency 
    % Y: frequency spectrum (2-sided, complex)
    L = size(y, 2);
    y = y';
    Y = fft(y);  
    w = Fs/L*(-L/2:L/2-1);
    P = abs(Y/L)'; P = P(:, 1:L/2+1);
    P(:, 2:end-1) = 2*P(:, 2:end-1);
    P = P.^2;
    wP = Fs/L*(0:(L/2));
    Y = fftshift(Y)';
end