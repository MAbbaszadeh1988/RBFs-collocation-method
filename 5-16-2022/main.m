clc
% h = waitbar(1,'Please wait...');
% k=0;
% tic
% for i=1:10000,
%     k=k+1;
%     waitbar(i/10000,h);
% end
% toc

  for i = 1:100000
      % Calculation
      workbar(i/100000,'Performing Calclations...','Progress') 
  end
