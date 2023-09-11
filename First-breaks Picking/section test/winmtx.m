function [ W ] = winmtx( N,sd )
%% -------->>> STFT Window Matrix Function
%   Detailed explanation goes here
W = zeros(N);
if length(sd) == 1
    W = toeplitz(fftshift(1/(sqrt(2*pi)*sd) * ...
        exp(-1/2*((-fix(N/2)+1:ceil(N/2))/sd).^2)));
else
    for i = 1:N
        W(:,i) = circshift(1/(sqrt(2*pi)*sd(i)) * ...
            exp(-1/2*((-fix(N/2):fix(N/2)-1)/sd(i)).^2)',fix(N/2)+i);
    end
end
end
