classdef Genetic_Al < handle
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        biType              % binary search only
        chromosomeNum       % total number of chromosomes
        totalGenration      % total generation
        stopCriteria        % stop criteria
        boundary            % boundary
        dimension           % dimension of the problem
        globalBest          % global best
        globalBestScore     % socre of global best
        globalBestRecord    % record all the value of global best
        totalResultScore    % total result
        ma                  % mather of offspring
        pa                  % father of offspring
        Chroms              % chromosomes 
        % target              % the target that we need to attend
        
        
    end
    
    methods
        function GA = Genetic_Al(totalGenration,dimension,boundary)
           % initilization of GA
           % for binary case the dimension will be structure dimensions
           % boundary will be the size of the structure
           GA.totalGenration = totalGenration;
           GA.dimension = dimension;           
           GA.boundary = boundary'; 
           if numel(boundary) == 2 && boundary(1) == 0 && boundary(2) == 1
                GA.biType = 1;        % specify optimization type
           end
           GA.globalBestScore = -inf;
           GA.globalBestRecord = zeros(totalGenration,3);
           %GA.globalBest = zeros(1,dimension);

        end
        
        function initialChromosome(GA,chromosomeNum,mutationRate,fitnessFnc)
            % initialization of chromosomes
            GA.chromosomeNum = chromosomeNum;
            for n = 1:chromosomeNum
                if ~isempty(GA.biType) 
                    pop = round(rand(1,GA.dimension)); % initilize population
                elseif isempty(GA.biType)
                    pop = GA.boundary(1,:) + (GA.boundary(2,:) - GA.boundary(1,:)).*rand(1,GA.dimension);
                else
                    error('check the simulation area')
                end
                bound = GA.boundary;
                if numel(mutationRate) == 1
                    C(n) = Chromosomes(n,pop,bound,mutationRate,fitnessFnc);
                elseif numel(mutationRate) == 2
                    randMutationNum = mutationRate(1) + (mutationRate(2) - mutationRate(1))*rand(1);
                    C(n) = Chromosomes(n,pop,bound,randMutationNum,fitnessFnc);
                else
                    error('check mutation input')
                end
            end
            GA.Chroms = C;
            GA.totalResultScore = zeros(1,chromosomeNum);
        end
        
        function evaluateFitFcn(GA,genNum)
            % this funcion evaluate all chromosomes
            totalValue = 0;
            for n = 1 : GA.chromosomeNum
                caculateFitness(GA.Chroms(n));
                GA.totalResultScore(n) = GA.Chroms(n).score;
                if GA.Chroms(n).score > GA.globalBestScore
                    GA.globalBestScore = GA.Chroms(n).score;
                    GA.globalBest = GA.Chroms(n).parameter;
                end
                GA.globalBestRecord(genNum,1) = genNum;
                GA.globalBestRecord(genNum,2) = GA.globalBestScore;
                totalValue = GA.Chroms(n).score + totalValue;
            end
            
            GA.globalBestRecord(genNum,3) = totalValue/GA.chromosomeNum;
        end
        
        function selectChrom(GA)
            % This function select and pass the best chromsome to next
            % generation (elitism method is used here)
            % find the best one
            [cost,ind]=sort(GA.totalResultScore,'descend'); % max element in first entry
            % GA.Chroms=GA.Chroms(ind,:); % sorts population with max cost first
            probs=cost/sum(cost); % simple normalization for probabilities
            matingNum=ceil(GA.chromosomeNum/2)-1; % number of matings
            mama=GA.RandChooseN(probs,matingNum); % mate #1
            papa=GA.RandChooseN(probs,matingNum); % mate #2
            GA.ma = ind(mama);
            GA.pa = ind(papa);
           
            
        end
        
        function Action=RandChooseN(GA,P,N)
            % function Action=RandChooseN(P,N)
            % Choose N numbers from 1 to length(P) using the
            % probabilities in P. For example, if P=[0.1,0.9],
            % we choose "1" 10% of the time, and "2" 90% of
            % the time. Selection is done WITH replacement,
            % so, for example, if N=3, we could return [2, 2, 2]
            %Set up the bins
            BinEdges=[0, cumsum(P(:)')];
            Action=zeros(1,N);
            for i=1:N
                x=rand;
                Counts=histc(x,BinEdges);
                Action(i)=find(Counts==1);
            end
        end
        
        function crossChrom(GA)
            % crossover and mutation of every chromosome
            matingNum=ceil(GA.chromosomeNum/2)-1; % number of matings
            offspring = ones(GA.chromosomeNum,prod(GA.dimension));
            for chromNum = 1 : matingNum
                father(chromNum,:) = GA.Chroms(GA.pa(chromNum)).parameter;
                mather(chromNum,:) = GA.Chroms(GA.ma(chromNum)).parameter;
                % single point corssover
%                 crossPoint = round(rand(1)*(prod(GA.dimension-1)));
%                 offspring(chromNum*2 -1,:) = [father(chromNum,1:crossPoint) mather(chromNum,crossPoint+1:end)];
%                 offspring(chromNum*2,:) = [mather(chromNum,1:crossPoint) father(chromNum,crossPoint+1:end)];
                % uniform crossover
                mask = round(rand(1,prod(GA.dimension)));
                offspring(chromNum*2 -1,:) = mask.*father(chromNum,:) + not(mask).*mather(chromNum,:);
                offspring(chromNum*2 ,:) = not(mask).*father(chromNum,:) + mask.*mather(chromNum,:);
            end
            offspring(matingNum*2+1,:) = GA.globalBest;
            offspring(matingNum*2+2,:) = GA.globalBest;
            for n = 1 : GA.chromosomeNum
                % give back the parameters for chromosomes
                GA.Chroms(n).setparameter(offspring(n,:))   
            end    
        end
        
        function mutateChrom(GA)
           % mutation of every chromosome
           for n = 1 : GA.chromosomeNum - 1
                GA.Chroms(n).mutateChromosome();  
            end    
        end
        
        
            
        function runGA(GA)
           % this is the main function of GA
           for geni = 1 : GA.totalGenration
               evaluateFitFcn(GA,geni);
               selectChrom(GA);
               crossChrom(GA);
               mutateChrom(GA);
           end
           
           
        end
        
        function PlotGbest(GA)
            % this function plot the result
            figure1 = figure;
            % Create axes
            axes1 = axes('Parent',figure1,'FontWeight','demi','FontSize',12);
            box(axes1,'on');
            hold(axes1,'all');
            plot(GA.globalBestRecord(:,1),GA.globalBestRecord(:,2),'--bo','linewidth',2);
            hold on;
            plot(GA.globalBestRecord(:,1),GA.globalBestRecord(:,3),'--r*','linewidth',2);
            hold off
            title('Evelution of the Global Best','FontWeight','bold','FontSize',14);
            xlabel('generations','FontWeight','demi','FontSize',12);
            ylabel('value of fitness function','FontWeight','demi','FontSize',12);
            legend('Global Best','Mean value');
        end
            
        
        
        
         
        
         
    end
    
end

