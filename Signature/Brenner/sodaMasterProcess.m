%% script sodaMasterProcess
% SODAMASTERPROCESS performs all processing steps and functions described
% in the processing document "sigProcessingGuide.pdf"
%
%   S.D.Brenner, 2019

% Clean workspace
clearvars;
clc;
close all hidden;

% List all moorings
sodaMoorings = {'A','B','C'};

% Loop through moorings and perform concatenation and processing steps
for s = 1:3
    % Get relevant Soda label
    sodaLabel = sodaMoorings{s};
    
    % Concatenate raw data
    % [ see processing document §4.1 ]
    sodaConcatRaw( sodaLabel );
    disp( ['SODA ',sodaLabel,' concatenation complete'] );
    
    % Initial QC Processing
    % [ see processing document §10 for summary, or §5-9 for details about
    %   specific steps ]
    sodaProcessSigAverage( sodaLabel );
    sodaProcessSigAverageIce( sodaLabel );
    sodaProcessSigBurst( sodaLabel );
    disp( ['SODA ',sodaLabel,' processing complete'] );
    
    % Build data products
    % [ These functions and products are still under-development ]
    sodaWavesProduct( sodaLabel );
    sodaIceProduct( sodaLabel );
    disp( ['SODA ',sodaLabel,' data products complete'] );
end