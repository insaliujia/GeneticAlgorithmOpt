% Genetic Algorithm
% The program is using object oriented program at MATLAB
% Jia LIU Ph.D student, INSA de Lyon
% Supervisor Regis Orobtchouk, INSA de Lyon
%% INITIALIZE MATLAB 
close all; 
% clear all; 
clc
% set(0,'ShowHiddenHandles','on'); delete(get(0,'Children'))
% addpath(genpath('D:\Gitcode\GA'));
addpath(genpath('D:\Gitcode\GeneticAlgorithmOpt'));
% define fitness function
Fitnessfnc = inline('sum(x,2)');
% define the simulation area 
% Genetic_Al(totalGenration,dimension,boundary)
GA = Genetic_Al(400,100,[0,1]);
% initialize chromosomes GA.initialChromosome(chromosomeNum,mutationRate,fitnessFnc)
GA.initialChromosome(40,0.01,Fitnessfnc);
% run GA
runGA(GA);
PlotGbest(GA)

