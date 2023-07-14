function [temp_lambda,temp_fare]=lower_model_final(delta,N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV)
max_lambda=12;
min_lambda=1;
%20201016max_fare=10;
max_fare=10;
%max_fare=2;
min_fare=0.1;
upper_lambda=max_lambda.*ones(1,M);
lower_lambda=min_lambda.*ones(1,M);
upper_fare=max_fare.*ones(1,M);
lower_fare=min_fare.*ones(1,M);
rng(3234);
initial_lambda=(max_lambda-min_lambda).*rand(1,M)+min_lambda;
initial_fare=(max_fare-min_fare).*rand(1,M)+min_fare;
% initial_lambda=min_lambda.*ones(1,M);
% initial_fare=min_fare.*ones(1,M);
sigma_sigma=1.5;
epsilon=0.001;
beta=0.5;
mu=2*ones(1,N+1);
dual_multi=ones(1,4*M+1);
[probability,demand,~,~]=cal_pro_demand_final(delta,initial_lambda,initial_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
[profit,~]=cal_profit(probability,demand,delta,M,N,common_true,distance_section_true,initial_fare,initial_lambda,max_distance,velocity,CV,dual_multi,min_lambda,max_lambda,min_fare,max_fare,mu,v);
%probability_temp=probability;
%demand_temp=demand;
temp_lambda=initial_lambda;
temp_fare=initial_fare;
temp_profit=profit;
%temp_L=L;
con=ones(1,4*M);
con(2:2:4*M)=-1;
constant=ones(1,4*M);
constant(1:2:2*M-1)=-min_lambda;
constant(2:2:2*M)=max_lambda;
constant(2*M+1:2:4*M-1)=-min_fare;
constant(2*M+2:2:4*M)=max_fare;
%aux_lambda=zeros(1,M);
%aux_fare=zeros(1,M);
k_iteration=1;
max_iteration=1000;
total_iteration=0;
while 1
    total_iteration=total_iteration+1;
    if total_iteration==51
        temp_lambda_upper=temp_lambda>upper_lambda;
        upper_lambda_index=find(temp_lambda_upper==true);
        if isempty(upper_lambda_index)==0
            temp_lambda(upper_lambda_index)=upper_lambda(upper_lambda_index);
        end
        temp_lambda_lower=temp_lambda<lower_lambda;
        lower_lambda_index=find(temp_lambda_lower==true);
        if isempty(lower_lambda_index)==0
            temp_lambda(lower_lambda_index)=lower_lambda(lower_lambda_index);
        end
        temp_fare_upper=temp_fare>upper_fare;
        upper_fare_index=find(temp_fare_upper==true);
        if isempty(upper_fare_index)==0
            temp_fare(upper_fare_index)=upper_fare(upper_fare_index);
        end
        temp_fare_lower=temp_fare<lower_fare;
        lower_fare_index=find(temp_fare_lower==true);
        if isempty(lower_fare_index)==0
            temp_fare(lower_fare_index)=lower_fare(lower_fare_index);
        end
        break;
    end
    %%
    lambda_temp=repmat(temp_lambda,2,1);
    lambda_temp=reshape(lambda_temp,1,[]);
    fare_temp=repmat(temp_fare,2,1);
    fare_temp=reshape(fare_temp,1,[]);
    variables=[lambda_temp,fare_temp];
    for t=1:max_iteration
        k_iteration=k_iteration+1;
        [probability_temp,demand_temp,~,~]=cal_pro_demand_final(delta,temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
        [partial_demand_lambda,partial_demand_fare,partial_probability_lambda]=partial(M,delta,temp_lambda,max_lambda,temp_fare,max_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,probability_temp,demand_temp,epsilon);
        [sum_partial_profit_lambda,sum_partial_profit_fare]=cal_partial_profit(probability_temp,demand_temp,partial_demand_lambda,partial_demand_fare,partial_probability_lambda,delta,M,N,temp_lambda,temp_fare,velocity,CV,v,common_true,distance_section_true,max_distance,dual_multi,mu,max_lambda,min_lambda,max_fare,min_fare);  
        aux_lambda=(sum_partial_profit_lambda<-epsilon).*upper_lambda+(sum_partial_profit_lambda>epsilon).*lower_lambda+(abs(sum_partial_profit_lambda)<=epsilon).*temp_lambda;
        aux_fare=(sum_partial_profit_fare<-epsilon).*upper_fare+(sum_partial_profit_fare>epsilon).*lower_fare+(abs(sum_partial_profit_fare)<=epsilon).*temp_fare;
        temp_lambda=temp_lambda+1/k_iteration.*(aux_lambda-temp_lambda);
        temp_fare=temp_fare+1/k_iteration.*(aux_fare-temp_fare);
    end
    %%
%     [partial_demand_lambda,partial_demand_fare,partial_probability_lambda]=partial(M,delta,temp_lambda,max_lambda,temp_fare,max_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,probability_temp,demand_temp);
%     [sum_partial_profit_lambda,sum_partial_profit_fare]=cal_partial_profit(probability_temp,demand_temp,partial_demand_lambda,partial_demand_fare,partial_probability_lambda,delta,M,N,temp_lambda,temp_fare,velocity,CV,v,potential_demand_true,common_true,distance_section_true,max_distance,dual_multi,mu,max_lambda,max_fare,sigma);
%     line_set=find(delta~=N+1);
%     if isempty(line_set)==0
%         aux_lambda(line_set)=(sum_partial_profit_lambda(line_set)<-epsilon).*lower_lambda(line_set)+(sum_partial_profit_lambda(line_set)>epsilon).*upper_lambda(line_set)+(abs(sum_partial_profit_lambda(line_set))<=epsilon).*sum_partial_profit_lambda(line_set);
%         aux_fare(line_set)=(sum_partial_profit_fare(line_set)<-epsilon).*lower_fare(line_set)+(sum_partial_profit_fare(line_set)>epsilon).*upper_fare(line_set)+(abs(sum_partial_profit_fare(line_set))<=epsilon).*sum_partial_profit_fare(line_set);
%     end
%     line_set_0=find(delta==N+1);
%     if isempty(line_set_0)==0
%         aux_lambda(line_set_0)=(sum_partial_profit_lambda(line_set_0)<-epsilon).*upper_lambda(line_set_0)+(sum_partial_profit_lambda(line_set_0)>epsilon).*lower_lambda(line_set_0)+(abs(sum_partial_profit_lambda(line_set_0))<=epsilon).*sum_partial_profit_lambda(line_set_0);
%         aux_fare(line_set_0)=(sum_partial_profit_fare(line_set_0)<-epsilon).*upper_fare(line_set_0)+(sum_partial_profit_fare(line_set_0)>epsilon).*lower_fare(line_set_0)+(abs(sum_partial_profit_fare(line_set_0))<=epsilon).*sum_partial_profit_fare(line_set_0);
%     end
%     temp_lambda=temp_lambda+1/k.*(aux_lambda-temp_lambda);
%     temp_fare=temp_fare+1/k.*(aux_fare-temp_fare);
   % [probability_temp,demand_temp,~,~]=cal_pro_demand_final(delta,temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
    %[partial_demand_lambda,partial_demand_fare,partial_probability_lambda]=partial(M,delta,temp_lambda,max_lambda,temp_fare,max_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,probability_temp,demand_temp);
    %[sum_partial_profit_lambda,sum_partial_profit_fare]=cal_partial_profit(probability_temp,demand_temp,partial_demand_lambda,partial_demand_fare,partial_probability_lambda,delta,M,N,temp_lambda,temp_fare,velocity,CV,v,potential_demand_true,common_true,distance_section_true,max_distance,dual_multi,mu,max_lambda,max_fare,sigma);
    %aux_lambda=(sum_partial_profit_lambda<-epsilon).*upper_lambda+(sum_partial_profit_lambda>epsilon).*lower_lambda+(abs(sum_partial_profit_lambda)<=epsilon).*temp_lambda;
    %aux_fare=(sum_partial_profit_fare<-epsilon).*upper_fare+(sum_partial_profit_fare>epsilon).*lower_fare+(abs(sum_partial_profit_fare)<=epsilon).*temp_fare;
%     temp_lambda=temp_lambda+1/k.*(aux_lambda-temp_lambda);
%     temp_fare=temp_fare+1/k.*(aux_fare-temp_fare);
%%
    old_lambda=temp_lambda;
    old_fare=temp_fare;
    [new_lambda,new_fare]=BFGS(old_lambda,max_lambda,min_lambda,old_fare,max_fare,min_fare,delta,M,N,epsilon,velocity,CV,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,dual_multi,mu);
    
%%
    %new_lambda,new_fare
    [probability_temp,demand_temp,~,~]=cal_pro_demand_final(delta,new_lambda,new_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
    [profit_1,~]=cal_profit(probability_temp,demand_temp,delta,M,N,common_true,distance_section_true,new_fare,new_lambda,max_distance,velocity,CV,dual_multi,min_lambda,max_lambda,min_fare,max_fare,mu,v);
    %criterion=sum(abs(profit_1(1:N)-temp_profit(1:N)),2)+abs(profit_1(N+1));
    %criterion=sum(abs(profit_1(1:N))-temp_profit(1:N),2)+sqrt(min(profit_1(N+1),0)^2);
    lambda_new=repmat(new_lambda,2,1);
    lambda_new=reshape(lambda_new,1,[]);
    fare_new=repmat(new_fare,2,1);
    fare_new=reshape(fare_new,1,[]);
    variables_new=[lambda_new,fare_new];
    criterion=sqrt(sum(min(con.*variables_new+constant,0).^2,2)+min(profit_1(N+1),0)^2);
    if criterion<=epsilon
        temp_lambda=new_lambda;
        temp_fare=new_fare;
        break;
    else
%         if isempty(line_set_0)==0
% %             if abs(profit_1(N+1))/abs(temp_profit(N+1))>beta
% %                 mu=sigma_sigma*mu;
% %             end
%             if sqrt(min(profit_1(N+1),0)^2)/sqrt(min(temp_profit(N+1),0)^2)>beta
%                 mu=sigma_sigma*mu;
%             end
%         end
        for n=1:N+1
            line_set=find(delta==n);
            index=zeros(1,4*size(line_set,2));
            if isempty(line_set)==0
                index_temp=repmat(line_set,2,1);
                index_temp(1,:)=2.*index_temp(1,:)-1;
                index_temp(2,:)=2.*index_temp(2,:);
                index_fare=2*M+index_temp;
                index(1:4*size(line_set,2))=[reshape(index_temp,1,[]),reshape(index_fare,1,[])];
                dual_multi(index)=max(dual_multi(index)-mu(n).*(con(index).*variables_new(index)+constant(index)),0);
                if n==N+1
                    dual_multi(end)=max(dual_multi(end)-mu(n)*profit_1(N+1),0);
                end
                if n==N+1
                    if sqrt(sum(min(con(index).*variables_new(index)+constant(index),0).^2,2)+min(profit_1(N+1),0)^2)/sqrt(sum(min(con(index).*variables(index)+constant(index),0).^2,2)+min(temp_profit(N+1),0)^2)>beta
                        mu(n)=sigma_sigma*mu(n);
                    end
                else
                    if sqrt(sum(min(con(index).*variables_new(index)+constant(index),0).^2,2))/sqrt(sum(min(con(index).*variables(index)+constant(index),0).^2,2))>beta
                        mu(n)=sigma_sigma*mu(n);
                    end
                end
            end            
        end
        %dual_multi=dual_multi-mu*profit_1(N+1);
%         dual_multi=max(dual_multi-mu*profit_1(N+1),0);
%         temp_profit=profit_1;
        temp_lambda=new_lambda;
        temp_fare=new_fare;
        temp_profit=profit_1;
    end
end
end