% Load data from Excel
dat = xlsread('Rajbhandari S-Project-W2425.xlsx');
NoTarget = xlsread('Rajbhandari S-Project-W2425.xlsx', 1, 'A1:A70'); % 70 samples for Target Absent
Target = xlsread('Rajbhandari S-Project-W2425.xlsx', 1, 'A71:A130'); % 60 samples for Target Present

%%Find Midpoint of Maximums
threshold_mid = ((xTarget(idm2))+xNoTarget(idm1))/2
%Use Midpoint as Threshold
CorrectNoTargetDet = NoTarget < threshold_mid;
NumCorrectNoTarget = sum(CorrectNoTargetDet);
CorrectTargetDet = Target > threshold_mid;
NumCorrectTarget = sum(CorrectTargetDet);
FalseAlarm = NoTarget > threshold_mid;
NumFalseAlarm = sum(FalseAlarm);
Miss = Target < threshold_mid;
NumMiss = sum(Miss);
%Confusion Matrix
ConfusionMatrix_mid = [NumCorrectNoTarget, NumMiss; NumFalseAlarm, NumCorrectTarget]
%Transition Matix
TransitionMatrix_mid = ConfusionMatrix_mid*[1/70,0;0,1/60]
%A priori probability
Apriori_mid=[70/130;60/130]
%Output Vector, product of transition matrix and apriori gives the output
%probabilities
OutputVector_mid = TransitionMatrix_mid*Apriori_mid
%Error Rate
ErrorRate_mid = (NumFalseAlarm + NumMiss)/130
%PPV
PPV_mid = NumCorrectTarget/(NumFalseAlarm+NumCorrectTarget)

%Use Intersection as Threshold
threshold_ = 3.9923
CorrectNoTargetDet = NoTarget < threshold_;
NumCorrectNoTarget = sum(CorrectNoTargetDet);
CorrectTargetDet = Target > threshold_;
NumCorrectTarget = sum(CorrectTargetDet);
FalseAlarm = NoTarget > threshold_;
NumFalseAlarm = sum(FalseAlarm);
Miss = Target < threshold_;
NumMiss = sum(Miss);
%Confusion Matrix
ConfusionMatrix_ = [NumCorrectNoTarget, NumMiss; NumFalseAlarm, NumCorrectTarget]
%Transition Matix
TransitionMatrix_ = ConfusionMatrix_*[1/70,0;0,1/60]
%A priori probability
Apriori_=[70/130;60/170]
%Output Vector, product of transition matrix and apriori gives the output
%probabilities
OutputVector_ = TransitionMatrix_*Apriori_
%Error Rate
ErrorRate_ = (NumFalseAlarm + NumMiss)/130
%PPV
PPV_ = NumCorrectTarget/(NumFalseAlarm+NumCorrectTarget)