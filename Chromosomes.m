classdef Chromosomes < handle
    %This class takes care of all the chromosomes
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        name        % name of chromosome
        score       % fitness of chromosome
        parameter   % parameter of the function
        geneticCode % genetic code of parameters
        boundary    % boundary of every parameter
        mutationRate;    % mutation rate
        fitnessFunction  % fitness function 
    end
    
    methods
        function chrom = Chromosomes(name,parameter,boundary,mutationRate,fitnessFnc)
        % initialization function 
            chrom.name = name;
            chrom.parameter = parameter;
            chrom.boundary = boundary;
            chrom.mutationRate = mutationRate;
            chrom.fitnessFunction = fitnessFnc;
        end
        
        function caculateFitness(chromosome)
            % this function caculate the fitness of chromosomes
            chromosome.score = chromosome.fitnessFunction(chromosome.parameter);
        end
        
        function setparameter(chromosome,para)
           % set chromosome parameter
           chromosome.parameter = para;
        end
        
        function mutateChromosome(chromosome)
           % This function mutate chromosome
           chromosomeChangeNum = ceil(numel(chromosome.parameter)*chromosome.mutationRate);
           chromosomeChangePoint = randi([1,numel(chromosome.parameter)],1,chromosomeChangeNum);
           chromosome.parameter(1,chromosomeChangePoint) = 1 - chromosome.parameter(1,chromosomeChangePoint); 
           
            
        end
        

        
    
    end
    
end

