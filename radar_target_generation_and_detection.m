clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
maxRange = 200;
rangeResolution = 1;
maxVel = 100;
rtt = 5.5;

%speed of light = 3e8
c = 3e8;
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
 
R = 110;
v = -20;

%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

B = c / (2 * rangeResolution);
Tchirp = rtt * 2 * maxRange / c; 
slope =  B / Tchirp;
                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %Tcransmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    % This is done inside the send_signal function
    
    % *%TODO* :
    %For each time sample we need update the Tcransmitted and
    %received signal. 
    [Tx(i), Rx(i)]  = signals(fc, slope, R, v, t(i));
    
    % *%TODO* :
    %Now by mixing the Tcransmit and Receive generate the beat signal
    %This is done by element wise maTcrix multiplication of Tcransmit and
    %Receiver Signal
    Mix(i) = Tx(i) .* Rx(i);
    
end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
mix_nd = reshape(Mix, [Nr, Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
range_fft = fft(mix_nd) ./ Nr;

 % *%TODO* :
% Take the absolute value of FFT output
range_fft = abs(range_fft);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the specTcrum.
% Hence we throw out half of the samples.
range_fft = range_fft(1:(Nr/2));

%plotting the range
figure ('Name','Range from First FFT')
subplot(3,1,1)

 % *%TODO* :
 % plot FFT output 

plot(range_fft); 
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = pow2db(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
subplot(3,1,2);
surf(doppler_axis,range_axis,RDM);
colorbar;
%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of training Cells in both the dimensions.
Tcr = 8;
Tcd = 2;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gcr = 4;
Gcd = 1;

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 9;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells

CUT = 1;
training_Gcrid = (2*Tcr+2*Gcr+1)*(2*Tcd+2*Gcd+1);
n_guard = (2*Gcr+1)*(2*Gcd+1) - CUT;
n_training = training_Gcrid - n_guard - CUT;

filtered_sig = zeros(size(RDM));

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.


   % Use RDM[x,y] as the maTcrix from the output of 2D FFT for implementing
   % CFAR

% we already know the noise leve
for i = Tcr + Gcr + 1 : Nr/2 - Tcr - Gcr
    for j = Tcd + Gcd + 1 : Nd - Tcd - Gcd
      train = db2pow(RDM(i - Tcr - Gcr : i + Tcr + Gcr, j - Tcd - Gcd : j + Tcd + Gcd));
      train(i - Gcr : i + Gcr, j - Gcd : j + Gcd) = 0;
      thresh = pow2db(sum(train, 'all') / n_training) + offset;
      if RDM(i,j) > thresh
        filtered_sig(i, j) = 1;
      end
    end
end

% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%maTcrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 
% Already done through initialization

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
subplot(3,1,3)
surf(doppler_axis,range_axis, filtered_sig);
colorbar;


 
 