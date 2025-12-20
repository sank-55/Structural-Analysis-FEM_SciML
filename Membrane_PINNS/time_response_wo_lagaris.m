
%% Without the Lagaris or Modified pde solution 

clc;
clear all;

% ------------------ DATA GENERATION ------------------

numBoundaryConditionPoints = [250 250 250 250];

x0BC1 = linspace(0,1,numBoundaryConditionPoints(1));
x0BC2 = linspace(0,1,numBoundaryConditionPoints(2));
x0BC3 = ones(1,numBoundaryConditionPoints(3));
x0BC4 = zeros(1,numBoundaryConditionPoints(4));

y0BC1 = zeros(1,numBoundaryConditionPoints(1));
y0BC2 = ones(1,numBoundaryConditionPoints(2));
y0BC3 = linspace(0,1,numBoundaryConditionPoints(3));
y0BC4 = linspace(0,1,numBoundaryConditionPoints(4));

c = 1;                              % wave speed
Tlimit = 1/(2*sqrt(2));

t0BC1 = linspace(0,Tlimit,numBoundaryConditionPoints(1));
t0BC2 = linspace(0,Tlimit,numBoundaryConditionPoints(2));
t0BC3 = linspace(0,Tlimit,numBoundaryConditionPoints(3));
t0BC4 = linspace(0,Tlimit,numBoundaryConditionPoints(4));

u0BC1 = zeros(1,numBoundaryConditionPoints(1));
u0BC2 = zeros(1,numBoundaryConditionPoints(2));
u0BC3 = zeros(1,numBoundaryConditionPoints(3));
u0BC4 = zeros(1,numBoundaryConditionPoints(4));

numInitialConditionPoints = 1000;

x0IC = linspace(0,1,numInitialConditionPoints);
y0IC = linspace(0,1,numInitialConditionPoints);
t0IC = zeros(1,numInitialConditionPoints);
u0IC = sin(pi*x0IC).*sin(pi*y0IC);

numInternalCollocationPoints = 5000;
points = lhsdesign(numInternalCollocationPoints,3);
dataT = Tlimit*points(:,3);
dataX = points(:,1);
dataY = points(:,2);

ds = arrayDatastore([dataX dataY dataT]);

% ------------------ NETWORK INITIALIZATION ------------------

numLayers  = 4;
numNeurons = 20;

parameters = struct;

sz = [numNeurons 3];
parameters.fc1.Weights = initializeHe(sz,3);
parameters.fc1.Bias    = initializeZeros([numNeurons 1]);  % 1st layer

for layerNumber = 2:numLayers-1
    name = "fc" + layerNumber;
    sz   = [numNeurons numNeurons];
    numIn = numNeurons;
    parameters.(name).Weights = initializeHe(sz,numIn);
    parameters.(name).Bias    = initializeZeros([numNeurons 1]);
end

sz = [1 numNeurons];
numIn = numNeurons;
parameters.("fc" + numLayers).Weights = initializeHe(sz,numIn);
parameters.("fc" + numLayers).Bias    = initializeZeros([1 1]);

% ------------------ TRAINING OPTIONS ------------------

numEpochs          = 5000;
miniBatchSize      = 1000;
executionEnvironment = "auto";
initialLearnRate   = 0.01;
decayRate          = 0.005;

% ------------------ DATALOADER ------------------

mbq = minibatchqueue(ds, ...
    'MiniBatchSize',miniBatchSize, ...
    'MiniBatchFormat','BC', ...
    'OutputEnvironment',executionEnvironment);

% IC and BC dlarrays (CB format: Channel x Batch)
dlX0IC = dlarray(x0IC,'CB');
dlY0IC = dlarray(y0IC,'CB');
dlT0IC = dlarray(t0IC,'CB');
dlU0IC = dlarray(u0IC,'CB');

dlX0BC1 = dlarray(x0BC1,'CB'); dlY0BC1 = dlarray(y0BC1,'CB'); dlT0BC1 = dlarray(t0BC1,'CB');
dlX0BC2 = dlarray(x0BC2,'CB'); dlY0BC2 = dlarray(y0BC2,'CB'); dlT0BC2 = dlarray(t0BC2,'CB');
dlX0BC3 = dlarray(x0BC3,'CB'); dlY0BC3 = dlarray(y0BC3,'CB'); dlT0BC3 = dlarray(t0BC3,'CB');
dlX0BC4 = dlarray(x0BC4,'CB'); dlY0BC4 = dlarray(y0BC4,'CB'); dlT0BC4 = dlarray(t0BC4,'CB');

dlU0BC1 = dlarray(u0BC1,'CB');
dlU0BC2 = dlarray(u0BC2,'CB');
dlU0BC3 = dlarray(u0BC3,'CB');
dlU0BC4 = dlarray(u0BC4,'CB');

% Optional: move IC/BC data to GPU
if (executionEnvironment == "auto" && canUseGPU) || (executionEnvironment == "gpu")
    dlX0IC = gpuArray(dlX0IC); dlY0IC = gpuArray(dlY0IC);
    dlT0IC = gpuArray(dlT0IC); dlU0IC = gpuArray(dlU0IC);
    dlX0BC1 = gpuArray(dlX0BC1); dlY0BC1 = gpuArray(dlY0BC1); dlT0BC1 = gpuArray(dlT0BC1);
    dlX0BC2 = gpuArray(dlX0BC2); dlY0BC2 = gpuArray(dlY0BC2); dlT0BC2 = gpuArray(dlT0BC2);
    dlX0BC3 = gpuArray(dlX0BC3); dlY0BC3 = gpuArray(dlY0BC3); dlT0BC3 = gpuArray(dlT0BC3);
    dlX0BC4 = gpuArray(dlX0BC4); dlY0BC4 = gpuArray(dlY0BC4); dlT0BC4 = gpuArray(dlT0BC4);
    dlU0BC1 = gpuArray(dlU0BC1); dlU0BC2 = gpuArray(dlU0BC2);
    dlU0BC3 = gpuArray(dlU0BC3); dlU0BC4 = gpuArray(dlU0BC4);
end

averageGrad   = [];
averageSqGrad = [];

figure(1)
C = colororder;
lineLoss = animatedline('Color',C(2,:));
ylim([0 inf])
xlabel("Iteration")
ylabel("Loss")
grid on

start = tic;
iteration = 0;

% ------------------ TRAINING LOOP ------------------

for epoch = 1:numEpochs
    reset(mbq);

    while hasdata(mbq)
        iteration = iteration + 1;

        dlXYT = next(mbq);              % 3 x batchSize, 'BC'
        dlX = dlXYT(1,:);               % 1 x B
        dlY = dlXYT(2,:);
        dlT = dlXYT(3,:);

        % Convert to 'CB' (channel x batch) for fullyconnect
        dlX = dlarray(dlX,'CB');
        dlY = dlarray(dlY,'CB');
        dlT = dlarray(dlT,'CB');

        if (executionEnvironment == "auto" && canUseGPU) || (executionEnvironment == "gpu")
            dlX = gpuArray(dlX); dlY = gpuArray(dlY); dlT = gpuArray(dlT);
        end

        [gradients,loss] = dlfeval(@modelGradients,parameters,dlX,dlY,dlT,...
            dlX0IC,dlY0IC,dlT0IC,dlU0IC,...
            dlX0BC1,dlY0BC1,dlT0BC1,...
            dlX0BC2,dlY0BC2,dlT0BC2,...
            dlX0BC3,dlY0BC3,dlT0BC3,...
            dlX0BC4,dlY0BC4,dlT0BC4,c);

        learningRate = initialLearnRate / (1+decayRate*iteration);

        [parameters,averageGrad,averageSqGrad] = adamupdate( ...
            parameters,gradients,averageGrad,averageSqGrad,iteration,learningRate);
    end

    lossPlot = double(gather(extractdata(loss)));
    addpoints(lineLoss,iteration, lossPlot);

    D = duration(0,0,toc(start),'Format','hh:mm:ss');
    title("Epoch: " + epoch + ", Elapsed: " + string(D) + ", Loss: " + lossPlot)
    drawnow
end

% ------------------ TESTING ------------------

tTest = [0.1 0.15 0.20 0.25];
numPredictions = 1000;
XTest_points = linspace(0,1,numPredictions);
YTest_points = linspace(0,1,numPredictions);
[Xmesh,Ymesh] = meshgrid(XTest_points,YTest_points);

for i = 1:length(tTest)
    t = tTest(i);
    TTest = t*ones(1,numPredictions);

    dlUPred = zeros(size(Xmesh),'like',Xmesh);

    for j = 1:size(Xmesh,1)
        XTest = Xmesh(j,:);
        YTest = Ymesh(j,:);

        dlXTest = dlarray(XTest,'CB');
        dlYTest = dlarray(YTest,'CB');
        dlTTest = dlarray(TTest,'CB');

        if (executionEnvironment == "auto" && canUseGPU) || (executionEnvironment == "gpu")
            dlXTest = gpuArray(dlXTest);
            dlYTest = gpuArray(dlYTest);
            dlTTest = gpuArray(dlTTest);
        end

        dlUPred(j,:) = modelU(parameters,dlXTest,dlYTest,dlTTest);
    end

    % UTest = sin(pi*Xmesh).*sin(pi*Ymesh)*cos(sqrt(2)*pi*t);
    % UTest_stor{i} = UTest;
    % 
    % UPred = extractdata(dlUPred);
    % UPred_stor{i} = UPred;
    % 
    % errU = UPred - UTest;
    % errU_stor{i} = errU;

     % exact solution
    UTest = sin(pi*Xmesh).*sin(pi*Ymesh)*cos(sqrt(2)*pi*t);
    UTest_stor{i} = UTest;

    % dlUPred is already a numeric array (double or gpuArray), no extractdata needed [web:43]
    UPred = dlUPred;
    UPred_stor{i} = UPred;

    errU = UPred - UTest;
    errU_stor{i} = errU;

    figure(2)
    subplot(2,2,i)
    surf(Xmesh,Ymesh,UPred,'FaceAlpha',0.5,'EdgeColor','none')
    zlim([0 1])
    colorbar
    title("Predicted response at t = " + t);

    figure(3)
    subplot(2,2,i)
    surf(Xmesh,Ymesh,UTest,'FaceAlpha',0.5,'EdgeColor','none')
    colorbar
    title("True response at t = " + t);

    figure(4)
    subplot(2,2,i)
    surf(Xmesh,Ymesh,errU,'FaceAlpha',0.5,'EdgeColor','none')
    colorbar
    title("Error at t = " + t);
end

% ------------------ LOSS / GRADIENT FUNCTION ------------------

function [gradients,loss] = modelGradients(parameters,dlX,dlY,dlT,...
    dlX0IC,dlY0IC,dlT0IC,dlU0IC,...
    dlX0BC1,dlY0BC1,dlT0BC1,...
    dlX0BC2,dlY0BC2,dlT0BC2,...
    dlX0BC3,dlY0BC3,dlT0BC3,...
    dlX0BC4,dlY0BC4,dlT0BC4,c)

    % Collocation predictions
    U = modelU(parameters,dlX,dlY,dlT);

    % First derivatives
    gradientsU = dlgradient(sum(U,'all'),{dlX,dlY,dlT}, ...
        'EnableHigherDerivatives',true);
    Ux = gradientsU{1};
    Uy = gradientsU{2};
    Ut = gradientsU{3};

    % Second derivatives
    Uxx = dlgradient(sum(Ux,'all'),dlX,'EnableHigherDerivatives',true);
    Uyy = dlgradient(sum(Uy,'all'),dlY,'EnableHigherDerivatives',true);
    Utt = dlgradient(sum(Ut,'all'),dlT,'EnableHigherDerivatives',true);

    % PDE residual: Utt - c^2 (Uxx + Uyy) = 0
    f1 = c*(Uxx + Uyy) - Utt;
    zeroTarget1 = zeros(size(f1), 'like', f1);
    lossF1 = mse(f1, zeroTarget1);
    lossF = lossF1;

    % Initial condition: u(x,y,0) = sin(pi x) sin(pi y)
    dlU0ICPred = modelU(parameters,dlX0IC,dlY0IC,dlT0IC);
    lossU0IC = mse(dlU0ICPred, dlU0IC);

    % Initial velocity: ut(x,y,0) = 0
    dlUdot0ICpred = dlgradient(sum(dlU0ICPred,'all'),dlT0IC, ...
        'EnableHigherDerivatives',true);
    Target9 = zeros(size(dlUdot0ICpred), 'like', dlUdot0ICpred);
    lossUdot0IC = mse(dlUdot0ICpred, Target9);

    % Boundary conditions: u = 0 on all boundaries
    dlU0BC1Pred = modelU(parameters,dlX0BC1,dlY0BC1,dlT0BC1);
    Target2 = zeros(size(dlU0BC1Pred),'like',dlU0BC1Pred);
    lossU0BC1 = mse(dlU0BC1Pred, Target2);

    dlU0BC2Pred = modelU(parameters,dlX0BC2,dlY0BC2,dlT0BC2);
    Target3 = zeros(size(dlU0BC2Pred),'like',dlU0BC2Pred);
    lossU0BC2 = mse(dlU0BC2Pred, Target3);

    dlU0BC3Pred = modelU(parameters,dlX0BC3,dlY0BC3,dlT0BC3);
    Target4 = zeros(size(dlU0BC3Pred),'like',dlU0BC3Pred);
    lossU0BC3 = mse(dlU0BC3Pred, Target4);

    dlU0BC4Pred = modelU(parameters,dlX0BC4,dlY0BC4,dlT0BC4);
    Target5 = zeros(size(dlU0BC4Pred),'like',dlU0BC4Pred);
    lossU0BC4 = mse(dlU0BC4Pred, Target5);

    lossU = lossU0IC + lossUdot0IC + lossU0BC1 + lossU0BC2 + lossU0BC3 + lossU0BC4;

    loss = lossF + lossU;

    gradients = dlgradient(loss,parameters);
end

% ------------------ MODEL FUNCTION ------------------

function dlU = modelU(parameters,dlX,dlY,dlT)
    dlXYT = [dlX; dlY; dlT];    % size: 3 x batch, 'CB'
    numLayers = numel(fieldnames(parameters));

    weights = parameters.fc1.Weights;
    bias    = parameters.fc1.Bias;
    dlU = fullyconnect(dlXYT,weights,bias);

    for i = 2:numLayers
        name = "fc" + i;
        dlU = sin(dlU);
        weights = parameters.(name).Weights;
        bias    = parameters.(name).Bias;
        dlU = fullyconnect(dlU, weights, bias);
    end
end
% He Initialization Function
function weights = initializeHe(sz, numIn)
    arguments
        sz
        numIn
    end

    % He initialization: weights ~ N(0, sqrt(2/numIn))
    if numel(sz) == 1
        sz = [sz, numIn];
    end

    stddev = sqrt(2/numIn);
    weights = randn(sz, 'single') * stddev;
    weights = dlarray(weights);
end

% Zeros Initialization Function
function parameter = initializeZeros(sz)
    parameter = zeros(sz, 'single');
    parameter = dlarray(parameter);
end


   
