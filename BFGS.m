function [new_lambda,new_fare]=BFGS(lambda,max_lambda,min_lambda,fare,max_fare,min_fare,delta,M,N,epsilon,velocity,CV,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,dual_multi,mu)
%probability,demand,partial_demand_lambda,partial_demand_fare,partial_probability_lambda,delta,M,N,final_h,final_f,velocity,CV,v,potential_demand_true,common_true,distance_section_true,max_distance,dual_multi,mu,max_h,min_h,max_f,min_f,sigma
index=zeros(1,M);
index_num=zeros(1,N+1);
H=cell(1,N+1);
for n=1:N+1
    line_set=find(delta==n);
    if isempty(line_set)==0
        line_set_size=size(line_set,2);
        index_num(1,n)=line_set_size;
        index(sum(index_num(1:n-1),2)+1:sum(index_num(1:n-1),2)+index_num(1,n))=line_set;
    end
    H{n}=eye(2*index_num(n));
end
k=ones(1,N+1);
stop_criterion=zeros(1,N+1);
% for n=1:N+1
%     H{n}=eye(index_num(n));
% end
old_lambda=lambda;
old_fare=fare;
total_iteration=1;
max_iteration=1000;
while 1
    if total_iteration==max_iteration+1
        new_lambda=old_lambda;
        new_fare=old_fare;
        break;
    end
    total_iteration=total_iteration+1;
    lambda_old=old_lambda;
    fare_old=old_fare;
    [probability,demand,~,~]=cal_pro_demand_final(delta,old_lambda,old_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
    [partial_demand_lambda,partial_demand_fare,partial_probability_lambda]=partial(M,delta,old_lambda,max_lambda,old_fare,max_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,probability,demand,epsilon);
    [sum_partial_profit_lambda,sum_partial_profit_fare]=cal_partial_profit(probability,demand,partial_demand_lambda,partial_demand_fare,partial_probability_lambda,delta,M,N,old_lambda,old_fare,velocity,CV,v,common_true,distance_section_true,max_distance,dual_multi,mu,max_lambda,min_lambda,max_fare,min_fare);
    [~,L]=cal_profit(probability,demand,delta,M,N,common_true,distance_section_true,old_fare,old_lambda,max_distance,velocity,CV,dual_multi,min_lambda,max_lambda,min_fare,max_fare,mu,v);
    gradient1=sum_partial_profit_lambda(index);
    gradient2=sum_partial_profit_fare(index);
    %lambda1=old_lambda(index);
    %fare1=old_fare(index);
    for n=1:N+1
        if index_num(n)~=0 && stop_criterion(n)==0
            g=[gradient1(sum(index_num(1:n-1),2)+1:sum(index_num(1:n-1),2)+index_num(1,n)),gradient2(sum(index_num(1:n-1),2)+1:sum(index_num(1:n-1),2)+index_num(1,n))];
            d=-H{n}*(g.');
            [new_alpha]=line_search(lambda_old,max_lambda,min_lambda,fare_old,max_fare,min_fare,delta,M,n,N,velocity,CV,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,dual_multi,mu,g,d,index,index_num);
            line_set=index(sum(index_num(1:n-1),2)+1:sum(index_num(1:n-1),2)+index_num(1,n));
            old_lambda(line_set)=old_lambda(line_set)+(new_alpha.*d(1:size(d,1)/2)).';
            old_fare(line_set)=old_fare(line_set)+(new_alpha.*d(size(d,1)/2+1:size(d,1))).';
        end
    end
    [probability_new,demand_new,~,~]=cal_pro_demand_final(delta,old_lambda,old_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
    [partial_demand_lambda_new,partial_demand_fare_new,partial_probability_lambda_new]=partial(M,delta,old_lambda,max_lambda,old_fare,max_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,probability_new,demand_new,epsilon);
    [sum_partial_profit_lambda_new,sum_partial_profit_fare_new]=cal_partial_profit(probability_new,demand_new,partial_demand_lambda_new,partial_demand_fare_new,partial_probability_lambda_new,delta,M,N,old_lambda,old_fare,velocity,CV,v,common_true,distance_section_true,max_distance,dual_multi,mu,max_lambda,min_lambda,max_fare,min_fare);
    [~,L_new]=cal_profit(probability_new,demand_new,delta,M,N,common_true,distance_section_true,old_fare,old_lambda,max_distance,velocity,CV,dual_multi,min_lambda,max_lambda,min_fare,max_fare,mu,v);
    for n=1:N+1
        if index_num(n)~=0
            line_set=index(sum(index_num(1:n-1),2)+1:sum(index_num(1:n-1),2)+index_num(1,n)); 
            if norm([sum_partial_profit_lambda_new(line_set),sum_partial_profit_fare_new(line_set)])<=epsilon || norm([old_lambda(line_set),old_fare(line_set)]-[lambda_old(line_set),fare_old(line_set)])<=1e-5 || norm(L_new(n)-L(n))<=1e-2
                stop_criterion(n)=1;
            end
        end
    end
    if sum(stop_criterion,2)==size(unique(delta),2)
        new_lambda=old_lambda;
        new_fare=old_fare;
        break;
    else
        for n=1:N+1
            if index_num(n)~=0 && stop_criterion(n)==0
                if k(n)==2*index_num(n)
                    H{n}=eye(2*index_num(n));
                    k(n)=0;
                else
                    line_set=index(sum(index_num(1:n-1),2)+1:sum(index_num(1:n-1),2)+index_num(1,n));
                    p_lambda=old_lambda(line_set)-lambda_old(line_set);
                    p_fare=old_fare(line_set)-fare_old(line_set);
                    q_lambda=sum_partial_profit_lambda_new(line_set)-sum_partial_profit_lambda(line_set);
                    q_fare=sum_partial_profit_fare_new(line_set)-sum_partial_profit_fare(line_set);
                    p=[p_lambda,p_fare];
                    q=[q_lambda,q_fare];
                    H{n}=H{n}+(1+q*H{n}*(q.')/(p*(q.')))*((p.')*p/(p*(q.')))-((p.')*q*H{n}+H{n}*(q.')*p)/(p*(q.'));
                end
            end
        end
    end
    index_count=find(index_num~=0);
    k(index_count)=k(index_count)+1;
end
end