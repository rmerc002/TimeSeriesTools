%% PARAMETERS
DATASET    = [];
UPPERBOUND = 791;
LOWERBOUND = 1;
STEPSIZE   = 1;
EXACT_PMP  = magicMP;

%% EXECUTION
Exact               = EXACT_PMP;
Exact(isnan(Exact)) = 0;

RMSE        = [];
I           = LOWERBOUND:STEPSIZE:UPPERBOUND;
% IDX         = I(BinarySplit(I));         % The order in which we should explore the subequence lengths
IDX         = I(MaxSplit(I,Exact));         % The order in which we should explore the subequence lengths
ApproxIDX   = zeros(size(Exact,1),1); % ApproxIDX(10) = 1 implies the 10th MP should be approximated by the 1st MP
f = waitbar(0.0,'Calculating ApproximatePMP');
for i = 1:length(IDX)
  waitbar(i/length(IDX),f,'Calculating ApproximatePMP');
  lengthIDX            = find(I == IDX(i)); % Get index of explored subsequence length
  ApproxIDX(lengthIDX) = lengthIDX;
  
  % Propogate update until we find an evaluated row
  for j = lengthIDX + 1 : length(ApproxIDX)
    if (j == ApproxIDX(j))
      break;
    end
    ApproxIDX(j) = lengthIDX;
  end
  Square = (Exact - Exact(ApproxIDX,:)) .^ 2;
  RMSE(end + 1) = sqrt(mean(Square(:)));
end
close(f);

%% GENERATE FIGURE
figure; hold on;
title('Rate of Approximation Convergence');
xticks(LOWERBOUND:STEPSIZE:UPPERBOUND);
xticks([0:0.10:1] * length(IDX));
xticklabels({'0', '0.10', '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90','1.0'});
yticks([0:0.25:1] * max(RMSE));
yticklabels({'0', '0.25', '0.50', '0.75', '1.0'});
ylim([0,max(RMSE)]);
xlabel('Fraction of Effort (Relative to Exact PMP Completion)');
ylabel('Accuracy norm(ApproxPMP - ExactPMP)');
plot(RMSE, 'Linewidth', 2);
hold off;