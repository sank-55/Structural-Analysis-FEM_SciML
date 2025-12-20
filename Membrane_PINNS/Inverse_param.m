%% INVERSE PROBLEM
%     Actual value of the Param in coefficient of the GDE is '1'
%

% clc;
% clear all;
% 
% % ------------------ DATA: COLLLOCATION POINTS ------------------
% Tlimit = 1/(2*sqrt(2));
% numInternalCollocationPoints = 25000;
% 
% points = lhsdesign(numInternalCollocationPoints,3);   % [x,y,t] in (0,1)^2 x (0,Tlimit)
% dataT = Tlimit*points(:,3);
% dataX = points(:,1);
% dataY = points(:,2);
% 
% ds = arrayDatastore([dataX dataY dataT]);
% 
% % ------------------ NETWORK DEFINITION ------------------
% numLayers  = 4;
% numNeurons = 20;
% 
% parameters = struct;
% 
% % first layer
% sz = [numNeurons 3];
% parameters.fc1.Weights = initializeHe(sz,3);
% parameters.fc1.Bias    = initializeZeros([numNeurons 1]);
% 
% % hidden layers
% for layerNumber = 2:numLayers-1
%     name = "fc" + layerNumber;
%     sz   = [numNeurons numNeurons];
%     numIn = numNeurons;
%     parameters.(name).Weights = initializeHe(sz,numIn);
%     parameters.(name).Bias    = initializeZeros([numNeurons 1]);
% end
% 
% % last layer
% sz = [1 numNeurons];
% numIn = numNeurons;
% lastName = "fc" + numLayers;
% parameters.(lastName).Weights = initializeHe(sz,numIn);
% parameters.(lastName).Bias    = initializeZeros([1 1]);
% 
% % identification parameter c (wave speed squared or coefficient)
% parameters.(lastName).opt_param = dlarray(1.25);   % learnable scalar [web:45][web:46]
% 
% % wrap weights and biases as dlarray
% paramNames = fieldnames(parameters);
% for k = 1:numel(paramNames)
%     layerName = paramNames{k};
%     parameters.(layerName).Weights = dlarray(parameters.(layerName).Weights);
%     parameters.(layerName).Bias    = dlarray(parameters.(layerName).Bias);
% end  % [web:9][web:43]
% 
% % ------------------ TRAINING OPTIONS ------------------
% numEpochs            = 500;  % NEED MAXIMUM NO. OF EPOCHS
% miniBatchSize        = 10000;
% executionEnvironment = "auto";
% initialLearnRate     = 0.01;
% decayRate            = 0.005;
% 
% % ------------------ MINIBATCHQUEUE ------------------
% mbq = minibatchqueue(ds, ...
%     'MiniBatchSize',miniBatchSize, ...
%     'MiniBatchFormat','BC', ...
%     'OutputEnvironment',executionEnvironment);  % [web:35][web:49]
% 
% averageGrad   = [];
% averageSqGrad = [];
% 
% figure(1)
% C = colororder;
% lineLoss = animatedline('Color',C(2,:));
% ylim([0 inf])
% xlabel("Iteration")
% ylabel("Loss")
% grid on
% 
% % simulated measurement data for inverse term
% n_data      = 5000;
% data_points = lhsdesign(n_data,3);
% t_data      = Tlimit*data_points(:,3);
% x_data      = data_points(:,1);
% y_data      = data_points(:,2);
% 
% figure(2)
% C2 = colororder;
% lineLoss2 = animatedline('Color',C2(2,:));
% ylim([0 inf])
% xlabel("Iteration")
% ylabel("Identified parameter")
% grid on
% 
% % move parameters to GPU if needed (after dlarray wrapping)
% if (executionEnvironment == "auto" && canUseGPU) || (executionEnvironment == "gpu")
%     for k = 1:numel(paramNames)
%         layerName = paramNames{k};
%         parameters.(layerName).Weights   = gpuArray(parameters.(layerName).Weights);
%         parameters.(layerName).Bias      = gpuArray(parameters.(layerName).Bias);
%         if isfield(parameters.(layerName),'opt_param')
%             parameters.(layerName).opt_param = gpuArray(parameters.(layerName).opt_param);
%         end
%     end
% end  % [web:6][web:43]
% 
% start = tic;
% iteration = 0;
% 
% loss_hist   = zeros(1,numEpochs);
% param_conv  = zeros(1,numEpochs);
% 
% % ------------------ TRAINING LOOP ------------------
% for epoch = 1:numEpochs
%     reset(mbq);
% 
%     while hasdata(mbq)
%         iteration = iteration + 1;
% 
%         dlXYT = next(mbq);      % size: 3 x batch, 'BC'
%         dlX = dlXYT(1,:);
%         dlY = dlXYT(2,:);
%         dlT = dlXYT(3,:);
% 
%         % convert to 'CB'
%         dlX = dlarray(dlX,'CB');
%         dlY = dlarray(dlY,'CB');
%         dlT = dlarray(dlT,'CB');
% 
%         if (executionEnvironment == "auto" && canUseGPU) || (executionEnvironment == "gpu")
%             dlX = gpuArray(dlX);
%             dlY = gpuArray(dlY);
%             dlT = gpuArray(dlT);
%         end
% 
%         % physics + inverse loss
%         [gradients,loss] = dlfeval(@modelGradients,parameters,numLayers,dlX,dlY,dlT, ...
%                                    t_data,x_data,y_data);  % [web:45][web:47]
% 
%         learningRate = initialLearnRate / (1+decayRate*iteration);
% 
%         [parameters,averageGrad,averageSqGrad] = adamupdate( ...
%             parameters,gradients,averageGrad,averageSqGrad,iteration,learningRate);
%     end
% 
%     lossPlot = double(gather(extractdata(loss)));
%     addpoints(lineLoss,iteration, lossPlot);
% 
%     D = duration(0,0,toc(start),'Format','hh:mm:ss');
%     title("Epoch: " + epoch + ", Elapsed: " + string(D) + ", Loss: " + lossPlot)
%     drawnow
% 
%     id_param = parameters.(lastName).opt_param;
%     loss2    = double(gather(extractdata(id_param)));
%     addpoints(lineLoss2,iteration, loss2);
% 
%     D2 = duration(0,0,toc(start),'Format','hh:mm:ss');
%     title("Epoch: " + epoch + ", Elapsed: " + string(D2) + ", Parameter: " + loss2)
%     drawnow
% 
%     loss_hist(epoch)  = lossPlot;
%     param_conv(epoch) = loss2;
% end
% 
% figure(4)
% plot(1:numEpochs,param_conv)
% xlabel('Epoch')
% ylabel('Identified parameter(True Value : 1')
% 
% % ------------------ MODEL GRADIENTS FUNCTION ------------------
% function [gradients,loss] = modelGradients(parameters,numLayers,dlX,dlY,dlT, ...
%                                            t_data,x_data,y_data)
% 
% % identification parameter
% c_update = parameters.("fc" + numLayers).opt_param;
% 
% % forward at collocation points
% U = modelU(parameters,dlX,dlY,dlT);
% 
% % first derivatives [web:9][web:19]
% gradientsU = dlgradient(sum(U,'all'),{dlX,dlY,dlT},'EnableHigherDerivatives',true);
% Ux = gradientsU{1};
% Uy = gradientsU{2};
% Ut = gradientsU{3};
% 
% % second derivatives
% Uxx = dlgradient(sum(Ux,'all'),dlX,'EnableHigherDerivatives',true);
% Uyy = dlgradient(sum(Uy,'all'),dlY,'EnableHigherDerivatives',true);
% Utt = dlgradient(sum(Ut,'all'),dlT,'EnableHigherDerivatives',true);
% 
% % PDE residual: c*(Uxx+Uyy) - Utt = 0
% f1 = c_update*(Uxx + Uyy) - Utt;
% zeroTarget1 = zeros(size(f1),'like',f1);
% lossF = mse(f1, zeroTarget1);
% 
% % ---------- inverse/data loss ----------
% t_data = t_data';
% x_data = x_data';
% y_data = y_data';
% 
% noise  = 0.02*rand(size(t_data));
% UTrue  = sin(pi*x_data).*sin(pi*y_data).*cos(sqrt(2)*pi*t_data) + noise;
% 
% dlTdata = dlarray(t_data,'CB');
% dlXdata = dlarray(x_data,'CB');
% dlYdata = dlarray(y_data,'CB');
% 
% if isa(dlX,'gpuArray')
%     dlTdata = gpuArray(dlTdata);
%     dlXdata = gpuArray(dlXdata);
%     dlYdata = gpuArray(dlYdata);
% end
% 
% UModel = modelU(parameters,dlXdata,dlYdata,dlTdata);
% dlUTrue = dlarray(UTrue,'CB');
% if isa(UModel,'gpuArray')
%     dlUTrue = gpuArray(dlUTrue);
% end
% 
% UDiffInv = UModel - dlUTrue;
% zeroTargetInv = zeros(size(UDiffInv),'like',UDiffInv);
% lossInv = mse(UDiffInv,zeroTargetInv);
% 
% % total loss
% loss = lossF + lossInv;
% 
% % gradients w.r.t all parameters including opt_param [web:45][web:9]
% gradients = dlgradient(loss,parameters);
% end
% 
% % ------------------ MODEL FUNCTION ------------------
% function dlU = modelU(parameters,dlX,dlY,dlT)
% dlXYT = [dlX; dlY; dlT];      % 3 x batch, 'CB'
% numLayersLocal = numel(fieldnames(parameters));
% 
% % first fully connected
% weights = parameters.fc1.Weights;
% bias    = parameters.fc1.Bias;
% dlU = fullyconnect(dlXYT,weights,bias);
% 
% % remaining layers
% for i = 2:numLayersLocal
%     name = "fc" + i;
%     if ~isfield(parameters.(name),'Weights')
%         continue;   % skip opt_param-only struct field
%     end
%     dlU = sin(dlU);
%     weights = parameters.(name).Weights;
%     bias    = parameters.(name).Bias;
%     dlU = fullyconnect(dlU, weights, bias);
% end
% end
% 
% % ------------------ INITIALIZATION HELPERS ------------------
% function W = initializeHe(sz,numIn)
% W = randn(sz,'single')*sqrt(2/numIn);
% end
% 
% function B = initializeZeros(sz)
% B = zeros(sz,'single');
% end
