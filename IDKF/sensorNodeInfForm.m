classdef sensorNodeInfForm < handle
    % Implementation of the IDKF, see Pfaff et al. "Information Form 
    % Distributed Kalman Filtering (IDKF) with Explicit Inputs"
    % @author Florian Pfaff
    % @date 2011-2021
   properties
      y=0;
      Ymat_glob=1;
      Cv=1;
      H=1;
      transformToEnd=1;
      operationsList={};
   end
   methods
      function obj = sensorNodeInfForm(x,P,R,H,P_est_list,transformToEnd)
        obj.Cv=R;
        obj.H=H;
        
        % Convert to information form
        obj.y=P\x;
        obj.Ymat_glob=sum(reshape(cell2mat(cellfun(@inv,P_est_list,'UniformOutput',false)),[size(P_est_list{1}),length(P_est_list)]),3);
        
        obj.operationsList=[{'Cinv_!lok!'},{'<trans>'},obj.operationsList];
        obj.transformToEnd=transformToEnd;
      end
      function predict(obj,u)
        arguments
            obj sensorNodeInfForm
            u (:,1) double = zeros(size(obj.y))
        end
        global A Cw %#ok<GVMIS> 
        Ymat_glob_predinv=A/obj.Ymat_glob*A'+Cw;
        obj.y=Ymat_glob_predinv\(A*(obj.Ymat_glob\obj.y)+u);
        
        obj.transformToEnd=obj.transformToEnd*obj.Ymat_glob/A; %used inversion properties to replace /(F/obj.Ymat_glob)
        
        obj.operationsList=[{'Ypred'},{'plus u'},{'F'},{'Yestinv'},{'<pred>'},obj.operationsList];
        
        obj.transformToEnd=obj.transformToEnd*Ymat_glob_predinv;
        obj.Ymat_glob=inv(Ymat_glob_predinv);
        if (max(abs(nonzeros(obj.transformToEnd)))>1e+10)&&all(abs(nonzeros(obj.transformToEnd))>1e+5) % This is to prevent a convergence to infinity
            obj.transformToEnd=obj.transformToEnd/4096;
        end
      end
      function filter(obj,z)
        global R_list H_list %#ok<GVMIS> 
        R_transformed_list=arrayfun(@(i){H_list{i}'/R_list{i}*H_list{i}},1:length(R_list));
        obj.y=obj.y+obj.H'/obj.Cv*z;
        obj.Ymat_glob=obj.Ymat_glob+sum(reshape(cell2mat(R_transformed_list),[size(R_transformed_list{1}),length(R_transformed_list)]),3);
        
        obj.operationsList=[{'plus z'},{'<filter>'},obj.operationsList];
      end
      function nodeRes=fuse(node1,node2)
          assert(isa(node2,'sensorNodeSetInfForm'));
          assert(isequal(node1.Ymat_glob,node2.Ymat_glob),'Ymat_glob should be equal. Maybe a different sequence or number of prediction or filter steps were performed?.');
          assert(isequal(node1.transformToEnd,node2.transformToEnd),'Incompatible transformToEnd.');
          nodeRes=sensorNodeSetInfForm(0,1,1,1,1,{1,1},1); % Init empty and overwrite
          nodeRes.y=node1.y+node2.y;
          nodeRes.Ymat_glob=node1.Ymat_glob;
          nodeRes.Yellipmat=minkExternalApprox(node1.Yellipmat,node2.Yellipmat,'trace',node1.transformToEnd);
          nodeRes.transformToEnd=node1.transformToEnd;
          nodeRes.operationsList=[{'<fuse>'},node1.operationsList];
          
          % Assume measurements model is as for node1
          nodeRes.R=node1.Cv;
          nodeRes.H=node1.H;
      end
      function state=x(obj)
        state=obj.Ymat_glob\obj.y;
      end
      function cov=C(obj)
        cov=inv(obj.Ymat_glob);
      end
   end
    
end 