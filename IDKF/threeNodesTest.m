% @author Florian Pfaff
% @date 2011-2021
% This file contains tests with 3 sensors
global A Cw Cv_list H_list %#ok<GVMIS> 
%% initialize
rng(0)
A=diag([1,2,3,4]);
tmp=rand(2);
Cv_list{1}=tmp'*tmp;
tmp=rand(2);
Cv_list{2}=tmp'*tmp;

H_list={[eye(2),zeros(2)],[zeros(2),eye(2)],[zeros(2,1),eye(2),zeros(2,1)]};
%     H_list={[eye(2),eye(2)],[eye(2),2*eye(2)]}; Works as well
[z1,z2]=deal(2*ones(2,1),5*ones(2,1));

x_s1=3*ones(4,1);
x_s2=ones(4,1);
x_s3=[1;2;3;4];
P_list={eye(4),[3,0.1,0,0;0.1,3,0,0;0,0,3,0.1;0,0,0.1,3],diag([1,2,3,4])};
Cw=eye(4);
u = [1;2;3;4];
%% 1 prediction + fusion
[x_central,P_central,K{2}]=fuse(x_s1,P_list{1},x_s2,P_list{2});
K{1}=eye(size(A))-K{2};
[x_central,P_central,K{4}]=fuse(x_central,P_central,x_s3,diag([1,2,3,4]));
K{3}=eye(size(A))-K{4};
%%
C_centralAfterFusion=P_central;
[x_central,P_central]=predictionStep(x_central,P_central,A,Cw,u);

transformToEnd=A*C_centralAfterFusion;

setSensor1=sensorNodeInfForm(x_s1,P_list{1},Cv_list{1},H_list{1},P_list,transformToEnd);
setSensor2=sensorNodeInfForm(x_s2,P_list{2},Cv_list{2},H_list{2},P_list,transformToEnd);
setSensor3=sensorNodeInfForm(x_s3,P_list{3},Cv_list{2},H_list{2},P_list,transformToEnd);

setSensor1.predict(u);
setSensor2.predict();
setSensor3.predict();

% Test result
x_tracksFused=setSensor1.Ymat_glob\(setSensor1.y+setSensor2.y+setSensor3.y);
C_tracksFused=inv(setSensor1.Ymat_glob);

assert(norm(x_central-x_tracksFused)<0.01);
assert(norm(P_central-C_tracksFused)<0.01);

    
function [x,C]=predictionStep(x,C,A,Cw,u)
    x=A*x+u;
    C=A*C*A'+Cw;
end

function [x,C,K]=updateStep(x,C,H,y,Cv)
    K=C*H'/(Cv+H*C*H');
	x=x+K*(y-H*x);
	C=C-K*H*C;
end

function [x,C,K]=fuse(x1,C1,x2,C2,varargin)
    [x,C,K]=updateStep(x1,C1,eye(size(C1)),x2,C2);
    for i=1:3:numel(varargin)
        [x,C,Kcurr]=updateStep(x,C,eye(size(C1)),x2,C2);
        K=cat(3,K,Kcurr);
    end
end
    