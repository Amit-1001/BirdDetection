function [ceps] = mfcc_common(input, fs, nceps)
% Function to extract MFCC features of a frame.
persistent melbank

input = input(:);	% Make the input a column vector
N = length(input);
fftsize = 512;
NofMelFilters = 29;
fs=44100; 

% if (nargin < 2); fs = 16000; end;

Ce1		= sum(input.^2);		% Frame energy
Ce		= log(max(Ce1,2e-22));	% floors to 2 X 10 raised to power -22

if isempty(melbank)
	melbank		= melbankm(NofMelFilters,fftsize,fs); % MelFilterBank is returned
end
%hamWindow	= 0.54 - 0.46*cos(2*pi*(0:N-1)/N);
     windowedinput = input;%.*hamWindow';
fftmag		= abs(fft(windowedinput,2*fftsize)).^2;
temp		= log10(melbank * fftmag(1:fftsize));
ceps 		= dct(temp);
ceps(1)		= Ce;
ceps		= ceps(1:nceps);
 
% To decorrelate MFC Coefficients
% We want a matrixdct(i,j) which is totalFilters x cepstralCoefficients
% in size. The i,j component is given by
%                cos( i * (j+0.5)/NofMelFilters pi )
% where we have assumed that i and j start at 0.

% NofCepCoeff = 13;
% mfccDCTMatrix = 1/sqrt(NofMelFilters/2)*cos((0:(NofCepCoeff-1))' * ...
%	(2*(0:(NofMelFilters-1))+1) * pi/2/NofMelFilters);
% mfccDCTMatrix(1,:) = mfccDCTMatrix(1,:) * sqrt(2)/2;
% ceps = mfccDCTMatrix * temp;
%%%%%% Above block is equivalent to the first NofCepCoeff coefficients of -> ceps = dct(temp);
return

%--------------------------------------------------------------------------

function [x]=melbankm(p,n,fs,fl,fh)
%MELBANKM determine matrix for a mel-spaced filterbank [X,MN,MX]=(P,N,FS,FL,FH)
%
% Inputs:   p   number of filters in filterbank
%           n   length of fft
%           fs  sample rate in Hz
%           fl  low end of the lowest filter as a fraction of fs (default = 0)
%           fh  high end of highest filter as a fraction of fs (default = 0.5)

% Outputs:  x   a sparse matrix containing the filterbank amplitudes
%               If x is the only output argument then size(x)=[p,1+floor(n/2)]
%               otherwise size(x)=[p,mx-mn+1]

% To plot filterbanks e.g.  plot(melbankm(20,256,8000)')

if nargin<5
	fh=0.5;
	if nargin<4
		fl=0;
	end
end
f0=700/fs;
fn2=floor(n/2);
lr=log((f0+fh)/(f0+fl))/(p+1);
% convert to fft bin numbers with 0 for DC term
bl=n*((f0+fl)*exp([0 1 p p+1]*lr)-f0);
b2=ceil(bl(2));
b3=floor(bl(3));
b1=floor(bl(1))+1;
b4=min(fn2,ceil(bl(4)))-1;
pf=log((f0+(b1:b4)/n)/(f0+fl))/lr;
fp=floor(pf);
pm=pf-fp;
k2=b2-b1+1;
k3=b3-b1+1;
k4=b4-b1+1;
r=[fp(k2:k4) 1+fp(1:k3)];
c=[k2:k4 1:k3];
v=2*[1-pm(k2:k4) pm(1:k3)];
x = sparse(r,c,v,p,n);
return
