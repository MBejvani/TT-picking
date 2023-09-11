function [ W ] = winmtx_1( N,sd )
%WINMTX Summary of this function goes here
%   Detailed explanation goes here
W = zeros(N);
if length(sd) == 1
    W = toeplitz(fftshift(1/(sqrt(2*pi)*sd) * ...
        exp(-1/2*((-fix(N/2):fix(N/2)-1)/sd).^2)));
else 
    for i = 1:N
        W(:,i) = circshift(1/(sqrt(2*pi)*sd(i)) * ...
            exp(-1/2*((-fix(N/2):fix(N/2)-1)/sd(i)).^2)',fix(N/2)+i);
    end
end


