function [P1,P2,f]=fftwrapper(play)
    % ------------------------------------------------------------------
    % Input             
    %   play    :   Boolean 
    %
    % Output
    %   P1      :   1D vector
    %   P2      :   1D vector
    %   f       :   1D vector
    %
    % Example usage
    %   [P1,P2,f]=fftwrapper(true)
    % ------------------------------------------------------------------
    
    N=80000;        % Number of samplings in sample1 and sample2
    shift1=190000;  % First index of sample1
    shift2=320000;  % First index of sample2
    fcutoff=700;    % Highest frequency in returned spectrum
    
    % Load sound file and convert from stereo to mono
    [ystereo,Fs]=wavread('lwrhuntrd-ns197.wav');
    deltat=1/Fs;
    ymono=(ystereo(:,1)+ystereo(:,2))/2;
    sample1=ymono(1+shift1:N+shift1);
    sample2=ymono(1+shift2:N+shift2);
    
    if (play==true)
        % Play sound file and samples
        disp('Playing full audio file...');
        sound(ystereo,Fs);
        disp('Playing sample 1...');
        sound(sample1,Fs);
        disp('Playing sample 2...');
        sound(sample2,Fs);
    end
    
    % Do FFTs
    p1=fft(sample1);
    p2=fft(sample2);
    P1=p1.*conj(p1);
    P2=p2.*conj(p2);
    f=linspace(0,N-1,N)./(N*deltat);
    
    % Crop vectors to the sizes we are interested in
    ifcutoff=find(abs(f-fcutoff)==min(abs(f-fcutoff)))-1;
    f=f(1:ifcutoff);
    P1=P1(1:ifcutoff);
    P2=P2(1:ifcutoff);
end