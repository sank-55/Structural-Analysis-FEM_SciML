%%
% Contains the Modified output of the Neural Network using the Boundary condition ( Lagaris Condition)
%

% clc; clear;
% 
% %% ================= SECTION 1: CONVENTIONAL PINN (FORWARD) ================
% fprintf('\n========== SECTION 1: CONVENTIONAL PINN - FORWARD SOLUTION ==========\n');
% 
% % --- Boundary condition points ---
% numBoundaryConditionPoints = [250 250 250 250];
% 
% x0BC1 = linspace(0,1,numBoundaryConditionPoints(1));  % y = 0
% x0BC2 = linspace(0,1,numBoundaryConditionPoints(2));  % y = 1
% x0BC3 = ones(1,numBoundaryConditionPoints(3));        % x = 1
% x0BC4 = zeros(1,numBoundaryConditionPoints(4));       % x = 0
% 
% y0BC1 = zeros(1,numBoundaryConditionPoints(1));
% y0BC2 = ones(1,numBoundaryConditionPoints(2));
% y0BC3 = linspace(0,1,numBoundaryConditionPoints(3));
% y0BC4 = linspace(0,1,numBoundaryConditionPoints(4));
% 
% Tlimit = 1/(2*sqrt(2));
% t0BC1 = linspace(0,Tlimit,numBoundaryConditionPoints(1));
% t0BC2 = linspace(0,Tlimit,numBoundaryConditionPoints(2));
% t0BC3 = linspace(0,Tlimit,numBoundaryConditionPoints(3));
% t0BC4 = linspace(0,Tlimit,numBoundaryConditionPoints(4));
% 
% u0BC1 = zeros(1,numBoundaryConditionPoints(1));
% u0BC2 = zeros(1,numBoundaryConditionPoints(2));
% u0BC3 = zeros(1,numBoundaryConditionPoints(3));
% u0BC4 = zeros(1,numBoundaryConditionPoints(4));
% 
% % --- Initial condition points ---
% numInitialConditionPoints = 1000;
% x0IC = linspace(0,1,numInitialConditionPoints);
% y0IC = linspace(0,1,numInitialConditionPoints);
% t0IC = zeros(1,numInitialConditionPoints);
% u0IC = sin(pi*x0IC).*sin(pi*y0IC);
% 
% % --- Internal collocation points ---
% numInternalCollocationPoints = 5000;
% points = lhsdesign(numInternalCollocationPoints,3);
% dataX = points(:,1);
% dataY = points(:,2);
% dataT = Tlimit*points(:,3);
% 
% ds = arrayDatastore([dataX,dataY,dataT],'IterationDimension',1);
% 
% % --- Network definition ---
% numLayers  = 4;
% numNeurons = 20;
% 
% parameters = struct;
% 
% % fc1: 3 inputs -> hidden
% sz = [numNeurons 3];
% parameters.fc1.Weights = initializeHe(sz,3);
% parameters.fc1.Bias    = initializeZeros([numNeurons 1]);
% 
% % hidden layers
% for layerNumber = 2:numLayers-1
%     name = ['fc',num2str(layerNumber)];
%     sz = [numNeurons numNeurons];
%     parameters.(name).Weights = initializeHe(sz,numNeurons);
%     parameters.(name).Bias    = initializeZeros([numNeurons 1]);
% end
% 
% % output layer: hidden -> 1
% sz = [1 numNeurons];
% parameters.(['fc',num2str(numLayers)]).Weights = initializeHe(sz,numNeurons);
% parameters.(['fc',num2str(numLayers)]).Bias    = initializeZeros([1 1]);
% 
% % --- Training options ---
% numEpochs        = 1000;
% miniBatchSize    = 1000;
% executionEnvironment = 'auto';
% initialLearnRate = 0.01;
% decayRate        = 0.005;
% 
% mbq = minibatchqueue(ds, ...
%     'MiniBatchSize',miniBatchSize, ...
%     'MiniBatchFormat','BC', ... % returns numeric, not formatted dlarray
%     'OutputEnvironment',executionEnvironment);
% 
% % IC/BC as unformatted dlarray (just numeric with tracing)
% dlX0IC = dlarray(x0IC(:));  % column
% dlY0IC = dlarray(y0IC(:));
% dlT0IC = dlarray(t0IC(:));
% dlU0IC = dlarray(u0IC(:));
% 
% dlX0BC1 = dlarray(x0BC1(:)); dlY0BC1 = dlarray(y0BC1(:)); dlT0BC1 = dlarray(t0BC1(:));
% dlX0BC2 = dlarray(x0BC2(:)); dlY0BC2 = dlarray(y0BC2(:)); dlT0BC2 = dlarray(t0BC2(:));
% dlX0BC3 = dlarray(x0BC3(:)); dlY0BC3 = dlarray(y0BC3(:)); dlT0BC3 = dlarray(t0BC3(:));
% dlX0BC4 = dlarray(x0BC4(:)); dlY0BC4 = dlarray(y0BC4(:)); dlT0BC4 = dlarray(t0BC4(:));
% 
% dlU0BC1 = dlarray(u0BC1(:));
% dlU0BC2 = dlarray(u0BC2(:));
% dlU0BC3 = dlarray(u0BC3(:));
% dlU0BC4 = dlarray(u0BC4(:));
% 
% if strcmp(executionEnvironment,'auto') && canUseGPU
%     executionEnvironment = 'gpu';
%     dlX0IC = gpuArray(dlX0IC); dlY0IC = gpuArray(dlY0IC);
%     dlT0IC = gpuArray(dlT0IC); dlU0IC = gpuArray(dlU0IC);
%     dlX0BC1 = gpuArray(dlX0BC1); dlY0BC1 = gpuArray(dlY0BC1); dlT0BC1 = gpuArray(dlT0BC1);
%     dlX0BC2 = gpuArray(dlX0BC2); dlY0BC2 = gpuArray(dlY0BC2); dlT0BC2 = gpuArray(dlT0BC2);
%     dlX0BC3 = gpuArray(dlX0BC3); dlY0BC3 = gpuArray(dlY0BC3); dlT0BC3 = gpuArray(dlT0BC3);
%     dlX0BC4 = gpuArray(dlX0BC4); dlY0BC4 = gpuArray(dlY0BC4); dlT0BC4 = gpuArray(dlT0BC4);
%     dlU0BC1 = gpuArray(dlU0BC1); dlU0BC2 = gpuArray(dlU0BC2);
%     dlU0BC3 = gpuArray(dlU0BC3); dlU0BC4 = gpuArray(dlU0BC4);
% end
% 
% averageGrad   = [];
% averageSqGrad = [];
% 
% figure('Name','Conventional PINN Training Loss');
% C = colororder;
% lineLoss = animatedline('Color',C(2,:));
% ylim([0 inf]); xlabel('Iteration'); ylabel('Loss'); grid on;
% 
% c = 1;
% start     = tic;
% iteration = 0;
% 
% for epoch = 1:numEpochs
%     reset(mbq);
%     while hasdata(mbq)
%         iteration = iteration + 1;
%         batch = next(mbq);        % numeric N x 3
%         x = batch(:,1);
%         y = batch(:,2);
%         t = batch(:,3);
% 
%         dlX = dlarray(x);        % N x 1
%         dlY = dlarray(y);
%         dlT = dlarray(t);
%         if strcmp(executionEnvironment,'gpu')
%             dlX = gpuArray(dlX); dlY = gpuArray(dlY); dlT = gpuArray(dlT);
%         end
% 
%         [gradients,loss] = dlfeval(@modelGradients_ConvPINN,parameters,dlX,dlY,dlT, ...
%             dlX0IC,dlY0IC,dlT0IC,dlU0IC, ...
%             dlX0BC1,dlY0BC1,dlT0BC1, ...
%             dlX0BC2,dlY0BC2,dlT0BC2, ...
%             dlX0BC3,dlY0BC3,dlT0BC3, ...
%             dlX0BC4,dlY0BC4,dlT0BC4,c);
% 
%         learningRate = initialLearnRate/(1 + decayRate*iteration);
%         [parameters,averageGrad,averageSqGrad] = adamupdate(parameters,gradients, ...
%             averageGrad,averageSqGrad,iteration,learningRate);
%     end
% 
%     lossVal = double(gather(extractdata(loss)));
%     addpoints(lineLoss,iteration,lossVal);
%     D = duration(0,0,toc(start),'Format','hh:mm:ss');
%     title(['Epoch ',num2str(epoch),', Elapsed: ',string(D),', Loss: ',num2str(lossVal)]);
%     drawnow;
% end
% 
% fprintf('Computing predictions for Conventional PINN...\n');
% 
% tTest          = [0.1 0.15 0.20 0.25];
% numPredictions = 200;
% XTestpoints    = linspace(0,1,numPredictions);
% YTestpoints    = linspace(0,1,numPredictions);
% [Xmesh,Ymesh]  = meshgrid(XTestpoints,YTestpoints);
% 
% UPredstori = cell(numel(tTest),1);
% UTeststori = cell(numel(tTest),1);
% errUstori  = cell(numel(tTest),1);
% 
% for i = 1:numel(tTest)
%     t = tTest(i);
%     x_flat = Xmesh(:);
%     y_flat = Ymesh(:);
%     t_flat = t*ones(size(x_flat));
% 
%     dlXTest = dlarray(x_flat);
%     dlYTest = dlarray(y_flat);
%     dlTTest = dlarray(t_flat);
%     if strcmp(executionEnvironment,'gpu')
%         dlXTest = gpuArray(dlXTest);
%         dlYTest = gpuArray(dlYTest);
%         dlTTest = gpuArray(dlTTest);
%     end
%     dlUPred = modelU_ConvPINN(parameters,dlXTest,dlYTest,dlTTest); % N x 1
%     UPred   = reshape(gather(extractdata(dlUPred)),size(Xmesh));
% 
%     UTest   = sin(pi*Xmesh).*sin(pi*Ymesh).*cos(sqrt(2)*pi*t);
%     errU    = UPred - UTest;
% 
%     UPredstori{i} = UPred;
%     UTeststori{i} = UTest;
%     errUstori{i}  = errU;
% end
% 
% figure('Name','Conventional PINN - Predicted Response');
% for i = 1:numel(tTest)
%     subplot(2,2,i);
%     surf(Xmesh,Ymesh,UPredstori{i},'FaceAlpha',0.5,'EdgeColor','none');
%     zlim([-1 1]); colorbar;
%     title(['Predicted response at t = ',num2str(tTest(i))]);
%     xlabel('x'); ylabel('y'); zlabel('u');
% end
% 
% figure('Name','Conventional PINN - Analytical Solution');
% for i = 1:numel(tTest)
%     subplot(2,2,i);
%     surf(Xmesh,Ymesh,UTeststori{i},'FaceAlpha',0.5,'EdgeColor','none');
%     zlim([-1 1]); colorbar;
%     title(['True response at t = ',num2str(tTest(i))]);
%     xlabel('x'); ylabel('y'); zlabel('u');
% end
% 
% figure('Name','Conventional PINN - Absolute Error');
% for i = 1:numel(tTest)
%     subplot(2,2,i);
%     surf(Xmesh,Ymesh,errUstori{i},'FaceAlpha',0.5,'EdgeColor','none');
%     colorbar;
%     title(['Error at t = ',num2str(tTest(i))]);
%     xlabel('x'); ylabel('y'); zlabel('error');
% end
% 
% fprintf('Conventional PINN complete.\n\n');
% 
% %% ================= SECTION 2: MODIFIED OUTPUT PINN (FORWARD) =============
% fprintf('========== SECTION 2: MODIFIED OUTPUT PINN - FORWARD SOLUTION ==========\n');
% clc;
% 
% Tlimit = 1/(2*sqrt(2));
% numInternalCollocationPoints = 5000;
% points = lhsdesign(numInternalCollocationPoints,3);
% dataX = points(:,1);
% dataY = points(:,2);
% dataT = Tlimit*points(:,3);
% ds    = arrayDatastore([dataX,dataY,dataT],'IterationDimension',1);
% 
% numLayers  = 4;
% numNeurons = 20;
% 
% parameters_mod = struct;
% sz = [numNeurons 3];
% parameters_mod.fc1.Weights = initializeHe(sz,3);
% parameters_mod.fc1.Bias    = initializeZeros([numNeurons 1]);
% for layerNumber = 2:numLayers-1
%     name = ['fc',num2str(layerNumber)];
%     sz = [numNeurons numNeurons];
%     parameters_mod.(name).Weights = initializeHe(sz,numNeurons);
%     parameters_mod.(name).Bias    = initializeZeros([numNeurons 1]);
% end
% sz = [1 numNeurons];
% parameters_mod.(['fc',num2str(numLayers)]).Weights = initializeHe(sz,numNeurons);
% parameters_mod.(['fc',num2str(numLayers)]).Bias    = initializeZeros([1 1]);
% 
% numEpochs        = 1000;
% miniBatchSize    = 1000;
% initialLearnRate = 0.01;
% decayRate        = 0.005;
% 
% mbq = minibatchqueue(ds, ...
%     'MiniBatchSize',miniBatchSize, ...
%     'MiniBatchFormat','BC', ...
%     'OutputEnvironment',executionEnvironment);
% 
% averageGrad   = [];
% averageSqGrad = [];
% 
% figure('Name','Modified PINN Training Loss');
% C = colororder;
% lineLoss = animatedline('Color',C(2,:));
% ylim([0 inf]); xlabel('Iteration'); ylabel('Loss'); grid on;
% 
% c = 1;
% start     = tic;
% iteration = 0;
% 
% for epoch = 1:numEpochs
%     reset(mbq);
%     while hasdata(mbq)
%         iteration = iteration + 1;
%         batch = next(mbq);
%         dlX   = dlarray(batch(:,1));
%         dlY   = dlarray(batch(:,2));
%         dlT   = dlarray(batch(:,3));
% 
%         [gradients,loss] = dlfeval(@modelGradients_ModifiedPINN,parameters_mod,dlX,dlY,dlT,c);
%         learningRate = initialLearnRate/(1 + decayRate*iteration);
%         [parameters_mod,averageGrad,averageSqGrad] = adamupdate(parameters_mod,gradients, ...
%             averageGrad,averageSqGrad,iteration,learningRate);
%     end
% 
%     lossVal = double(gather(extractdata(loss)));
%     addpoints(lineLoss,iteration,lossVal);
%     D = duration(0,0,toc(start),'Format','hh:mm:ss');
%     title(['Epoch ',num2str(epoch),', Elapsed: ',string(D),', Loss: ',num2str(lossVal)]);
%     drawnow;
% end
% 
% fprintf('Computing predictions for Modified Output PINN...\n');
% 
% tTest          = [0.1 0.15 0.20 0.25];
% numPredictions = 200;
% XTestpoints    = linspace(0,1,numPredictions);
% YTestpoints    = linspace(0,1,numPredictions);
% [Xmesh,Ymesh]  = meshgrid(XTestpoints,YTestpoints);
% 
% UPredstori = cell(numel(tTest),1);
% UTeststori = cell(numel(tTest),1);
% errUstori  = cell(numel(tTest),1);
% 
% for i = 1:numel(tTest)
%     t = tTest(i);
%     x_flat = Xmesh(:);
%     y_flat = Ymesh(:);
%     t_flat = t*ones(size(x_flat));
% 
%     dlXTest = dlarray(x_flat);
%     dlYTest = dlarray(y_flat);
%     dlTTest = dlarray(t_flat);
%     if strcmp(executionEnvironment,'gpu')
%         dlXTest = gpuArray(dlXTest);
%         dlYTest = gpuArray(dlYTest);
%         dlTTest = gpuArray(dlTTest);
%     end
%     dlUPred = modelU_ModifiedPINN(parameters_mod,dlXTest,dlYTest,dlTTest);
%     UPred   = reshape(gather(extractdata(dlUPred)),size(Xmesh));
% 
%     UTest   = sin(pi*Xmesh).*sin(pi*Ymesh).*cos(sqrt(2)*pi*t);
%     errU    = UPred - UTest;
% 
%     UPredstori{i} = UPred;
%     UTeststori{i} = UTest;
%     errUstori{i}  = errU;
% end
% 
% figure('Name','Modified PINN - Predicted Response');
% for i = 1:numel(tTest)
%     subplot(2,2,i);
%     surf(Xmesh,Ymesh,UPredstori{i},'FaceAlpha',0.5,'EdgeColor','none');
%     zlim([-1 1]); colorbar;
%     title(['Predicted response at t = ',num2str(tTest(i))]);
%     xlabel('x'); ylabel('y'); zlabel('u');
% end
% 
% figure('Name','Modified PINN - Analytical Solution');
% for i = 1:numel(tTest)
%     subplot(2,2,i);
%     surf(Xmesh,Ymesh,UTeststori{i},'FaceAlpha',0.5,'EdgeColor','none');
%     zlim([-1 1]); colorbar;
%     title(['True response at t = ',num2str(tTest(i))]);
%     xlabel('x'); ylabel('y'); zlabel('u');
% end
% 
% figure('Name','Modified PINN - Absolute Error');
% for i = 1:numel(tTest)
%     subplot(2,2,i);
%     surf(Xmesh,Ymesh,errUstori{i},'FaceAlpha',0.5,'EdgeColor','none');
%     colorbar;
%     title(['Error at t = ',num2str(tTest(i))]);
%     xlabel('x'); ylabel('y'); zlabel('error');
% end
% 
% fprintf('Modified Output PINN complete.\n\n');
% 
% %% ================= SECTION 3: INVERSE PINN (IDENTIFY c) ==================
% fprintf('========== SECTION 3: INVERSE PINN - PARAMETER IDENTIFICATION ==========\n');
% clc;
% 
% Tlimit = 1/(2*sqrt(2));
% 
% numInternalCollocationPoints = 25000;
% points = lhsdesign(numInternalCollocationPoints,3);
% dataX = points(:,1);
% dataY = points(:,2);
% dataT = Tlimit*points(:,3);
% ds    = arrayDatastore([dataX,dataY,dataT],'IterationDimension',1);
% 
% numLayers  = 4;
% numNeurons = 20;
% 
% parameters_inv = struct;
% sz = [numNeurons 3];
% parameters_inv.fc1.Weights = initializeHe(sz,3);
% parameters_inv.fc1.Bias    = initializeZeros([numNeurons 1]);
% for layerNumber = 2:numLayers-1
%     name = ['fc',num2str(layerNumber)];
%     sz = [numNeurons numNeurons];
%     parameters_inv.(name).Weights = initializeHe(sz,numNeurons);
%     parameters_inv.(name).Bias    = initializeZeros([numNeurons 1]);
% end
% sz = [1 numNeurons];
% parameters_inv.(['fc',num2str(numLayers)]).Weights = initializeHe(sz,numNeurons);
% parameters_inv.(['fc',num2str(numLayers)]).Bias    = initializeZeros([1 1]);
% 
% % trainable wave speed c
% parameters_inv.(['fc',num2str(numLayers)]).optparam = dlarray(1.25);
% 
% numEpochs        = 3000;
% miniBatchSize    = 10000;
% initialLearnRate = 0.01;
% decayRate        = 0.005;
% 
% mbq = minibatchqueue(ds, ...
%     'MiniBatchSize',miniBatchSize, ...
%     'MiniBatchFormat','BC', ...
%     'OutputEnvironment',executionEnvironment);
% 
% averageGrad   = [];
% averageSqGrad = [];
% 
% figure('Name','Inverse PINN - Training Loss');
% C = colororder;
% lineLoss  = animatedline('Color',C(2,:));
% ylim([0 inf]); xlabel('Iteration'); ylabel('Total Loss'); grid on;
% 
% figure('Name','Inverse PINN - Identified Parameter c');
% C2 = colororder;
% lineLoss2 = animatedline('Color',C2(2,:));
% ylim([0 2]); xlabel('Iteration'); ylabel('Identified c'); grid on;
% 
% c_true = 1;
% ndata  = 5000;
% datapoints = lhsdesign(ndata,3);
% tdata = Tlimit*datapoints(:,1);
% xdata = datapoints(:,2);
% ydata = datapoints(:,3);
% 
% noise_level = 0.02;
% UTrue = sin(pi*xdata).*sin(pi*ydata).*cos(sqrt(2)*pi*tdata) + ...
%         noise_level*(rand(size(tdata))-0.5);
% 
% dlTdata = dlarray(tdata);
% dlXdata = dlarray(xdata);
% dlYdata = dlarray(ydata);
% dlUTrue = dlarray(UTrue);
% if strcmp(executionEnvironment,'gpu')
%     dlTdata = gpuArray(dlTdata);
%     dlXdata = gpuArray(dlXdata);
%     dlYdata = gpuArray(dlYdata);
%     dlUTrue = gpuArray(dlUTrue);
% end
% 
% start     = tic;
% iteration = 0;
% 
% for epoch = 1:numEpochs
%     reset(mbq);
%     while hasdata(mbq)
%         iteration = iteration + 1;
%         batch = next(mbq);
%         dlX   = dlarray(batch(:,1));
%         dlY   = dlarray(batch(:,2));
%         dlT   = dlarray(batch(:,3));
% 
%         [gradients,loss] = dlfeval(@modelGradients_InversePINN,parameters_inv,numLayers, ...
%             dlX,dlY,dlT,dlTdata,dlXdata,dlYdata,dlUTrue);
% 
%         learningRate = initialLearnRate/(1 + decayRate*iteration);
%         [parameters_inv,averageGrad,averageSqGrad] = adamupdate(parameters_inv,gradients, ...
%             averageGrad,averageSqGrad,iteration,learningRate);
%     end
% 
%     lossVal = double(gather(extractdata(loss)));
%     addpoints(lineLoss,iteration,lossVal);
%     D = duration(0,0,toc(start),'Format','hh:mm:ss');
%     title(['Epoch ',num2str(epoch),', Elapsed: ',string(D),', Loss: ',num2str(lossVal)]);
%     drawnow;
% 
%     c_identified = parameters_inv.(['fc',num2str(numLayers)]).optparam;
%     c_val        = double(gather(extractdata(c_identified)));
%     addpoints(lineLoss2,iteration,c_val);
%     drawnow;
% end
% 
% fprintf('\nInverse Problem Results:\n');
% c_final = double(gather(extractdata(parameters_inv.(['fc',num2str(numLayers)]).optparam)));
% fprintf('True wave speed c:       %.6f\n',c_true);
% fprintf('Identified wave speed c: %.6f\n',c_final);
% fprintf('Error in identification: %.6f\n\n',abs(c_final-c_true));
% 
% fprintf('All Example 4.4 computations complete!\n');
% 
% %% =========================== HELPER FUNCTIONS ============================
% 
% function W = initializeHe(sz,numIn)
%     bound = sqrt(2/numIn);
%     W = (2*bound)*rand(sz) - bound;
% end
% 
% function Z = initializeZeros(sz)
%     Z = zeros(sz);
% end
% 
% %% MODEL: Conventional PINN  (X = [dlX dlY dlT], all unformatted)
% function dlU = modelU_ConvPINN(parameters,dlX,dlY,dlT)
%     X = [dlX dlY dlT];            % N x 3
%     numLayers = numel(fieldnames(parameters));
% 
%     W = parameters.fc1.Weights;
%     B = parameters.fc1.Bias;
%     dlU = fullyconnect(X,W,B);    % N x numNeurons
% 
%     for i = 2:numLayers
%         dlU = tanh(dlU);
%         name = ['fc',num2str(i)];
%         W = parameters.(name).Weights;
%         B = parameters.(name).Bias;
%         dlU = fullyconnect(dlU,W,B); % last: N x 1
%     end
% end
% 
% %% MODEL: Modified output PINN
% function dlU = modelU_ModifiedPINN(parameters,dlX,dlY,dlT)
%     X = [dlX dlY dlT];  % N x 3
%     numLayers = numel(fieldnames(parameters));
% 
%     W = parameters.fc1.Weights;
%     B = parameters.fc1.Bias;
%     net = fullyconnect(X,W,B);
% 
%     for i = 2:numLayers
%         net = tanh(net);
%         name = ['fc',num2str(i)];
%         W = parameters.(name).Weights;
%         B = parameters.(name).Bias;
%         net = fullyconnect(net,W,B);
%     end
% 
%     x = dlX; y = dlY; t = dlT;
%     dlU = net .* (t.^2) .* (x.*(1-x)) .* (y.*(1-y)) .* sin(pi*x).*sin(pi*y);
% end
% 
% %% MODEL: Inverse PINN (network only; c is separate parameter)
% function dlU = modelU_InversePINN(parameters,dlX,dlY,dlT)
%     X = [dlX dlY dlT];  % N x 3
% 
%     fn = fieldnames(parameters);
%     numLayers = 0;
%     for k = 1:numel(fn)
%         if startsWith(fn{k},'fc') && isfield(parameters.(fn{k}),'Weights')
%             numLayers = numLayers + 1;
%         end
%     end
% 
%     W = parameters.fc1.Weights;
%     B = parameters.fc1.Bias;
%     dlU = fullyconnect(X,W,B);
% 
%     for i = 2:numLayers
%         dlU = tanh(dlU);
%         name = ['fc',num2str(i)];
%         W = parameters.(name).Weights;
%         B = parameters.(name).Bias;
%         dlU = fullyconnect(dlU,W,B);
%     end
% end
% 
% %% GRADIENT: Conventional PINN
% function [gradients,loss] = modelGradients_ConvPINN(parameters,dlX,dlY,dlT, ...
%     dlX0IC,dlY0IC,dlT0IC,dlU0IC, ...
%     dlX0BC1,dlY0BC1,dlT0BC1, ...
%     dlX0BC2,dlY0BC2,dlT0BC2, ...
%     dlX0BC3,dlY0BC3,dlT0BC3, ...
%     dlX0BC4,dlY0BC4,dlT0BC4,c)
% 
%     U = modelU_ConvPINN(parameters,dlX,dlY,dlT);   % N x 1
% 
%     g = dlgradient(sum(U,'all'),{dlX,dlY,dlT},'EnableHigherDerivatives',true);
%     Ux = g{1};
%     Uy = g{2};
%     Ut = g{3};
% 
%     Uxx = dlgradient(sum(Ux,'all'),dlX,'EnableHigherDerivatives',true);
%     Uyy = dlgradient(sum(Uy,'all'),dlY,'EnableHigherDerivatives',true);
%     Utt = dlgradient(sum(Ut,'all'),dlT,'EnableHigherDerivatives',true);
% 
%     f1    = c^2*(Uxx + Uyy) - Utt;
%     lossF = mse(f1,zeros(size(f1),'like',f1));
% 
%     % IC
%     U0Pred   = modelU_ConvPINN(parameters,dlX0IC,dlY0IC,dlT0IC);
%     lossU0IC = mse(U0Pred,dlU0IC);
% 
%     U0dot    = dlgradient(sum(U0Pred,'all'),dlT0IC,'EnableHigherDerivatives',true);
%     lossUdot = mse(U0dot,zeros(size(U0dot),'like',U0dot));
% 
%     % BCs
%     U0BC1 = modelU_ConvPINN(parameters,dlX0BC1,dlY0BC1,dlT0BC1);
%     U0BC2 = modelU_ConvPINN(parameters,dlX0BC2,dlY0BC2,dlT0BC2);
%     U0BC3 = modelU_ConvPINN(parameters,dlX0BC3,dlY0BC3,dlT0BC3);
%     U0BC4 = modelU_ConvPINN(parameters,dlX0BC4,dlY0BC4,dlT0BC4);
% 
%     lossBC1 = mse(U0BC1,zeros(size(U0BC1),'like',U0BC1));
%     lossBC2 = mse(U0BC2,zeros(size(U0BC2),'like',U0BC2));
%     lossBC3 = mse(U0BC3,zeros(size(U0BC3),'like',U0BC3));
%     lossBC4 = mse(U0BC4,zeros(size(U0BC4),'like',U0BC4));
% 
%     lossU = lossU0IC + lossUdot + lossBC1 + lossBC2 + lossBC3 + lossBC4;
% 
%     loss      = lossF + lossU;
%     gradients = dlgradient(loss,parameters);
% end
% 
% %% GRADIENT: Modified PINN
% function [gradients,loss] = modelGradients_ModifiedPINN(parameters,dlX,dlY,dlT,c)
%     U = modelU_ModifiedPINN(parameters,dlX,dlY,dlT);
% 
%     g = dlgradient(sum(U,'all'),{dlX,dlY,dlT},'EnableHigherDerivatives',true);
%     Ux = g{1};
%     Uy = g{2};
%     Ut = g{3};
% 
%     Uxx = dlgradient(sum(Ux,'all'),dlX,'EnableHigherDerivatives',true);
%     Uyy = dlgradient(sum(Uy,'all'),dlY,'EnableHigherDerivatives',true);
%     Utt = dlgradient(sum(Ut,'all'),dlT,'EnableHigherDerivatives',true);
% 
%     f1  = c^2*(Uxx + Uyy) - Utt;
%     loss = mse(f1,zeros(size(f1),'like',f1));
% 
%     gradients = dlgradient(loss,parameters);
% end
% 
% %% GRADIENT: Inverse PINN
% function [gradients,loss] = modelGradients_InversePINN(parameters,numLayers, ...
%     dlX,dlY,dlT,dlTdata,dlXdata,dlYdata,dlUTrue)
% 
%     c_update = parameters.(['fc',num2str(numLayers)]).optparam;
% 
%     U = modelU_InversePINN(parameters,dlX,dlY,dlT);
% 
%     g = dlgradient(sum(U,'all'),{dlX,dlY,dlT},'EnableHigherDerivatives',true);
%     Ux = g{1};
%     Uy = g{2};
%     Ut = g{3};
% 
%     Uxx = dlgradient(sum(Ux,'all'),dlX,'EnableHigherDerivatives',true);
%     Uyy = dlgradient(sum(Uy,'all'),dlY,'EnableHigherDerivatives',true);
%     Utt = dlgradient(sum(Ut,'all'),dlT,'EnableHigherDerivatives',true);
% 
%     f1    = c_update^2*(Uxx + Uyy) - Utt;
%     lossF = mse(f1,zeros(size(f1),'like',f1));
% 
%     UModel = modelU_InversePINN(parameters,dlXdata,dlYdata,dlTdata);
%     UDiff  = UModel - dlUTrue;
%     lossInv = mse(UDiff,zeros(size(UDiff),'like',UDiff));
% 
%     loss      = lossF + lossInv;
%     gradients = dlgradient(loss,parameters);
% end
