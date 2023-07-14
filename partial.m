 function [partial_demand_lambda,partial_demand_fare,partial_probability_lambda]=partial(M,delta,final_h,max_h,final_f,max_f,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,probability,demand,epsilon)
%function [partial_demand_lambda,partial_demand_fare,partial_probability_lambda]=partial(M,delta,final_h,final_f,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,probability,demand)
partial_demand_lambda=zeros(size(common_true,1),M);
partial_demand_fare=zeros(size(common_true,1),M);
partial_probability_lambda=zeros(size(common_true,1),M,M);
error=0.00001;
%epsilon=0.0001;
%min_lambda=min_h*ones(1,M);
max_lambda=max_h*ones(1,M);
%min_fare=min_f*ones(1,M);
max_fare=max_f*ones(1,M);
%[probability,demand,EW,ET]=cal_pro_demand_final(delta,final_h,final_f,velocity,sigma,v,potential_demand_true,common_true,distance_section_true)
for m=1:M
    temp_h=final_h;
    temp_f=final_f;
    if abs(final_h(1,m)-max_lambda(1,m))<=epsilon
        temp_h(1,m)=final_h(1,m)-error;
    else
        temp_h(1,m)=final_h(1,m)+error;
    end
    if abs(final_f(1,m)-max_fare(1,m))<=epsilon
        temp_f(1,m)=final_f(1,m)-error;
    else
        temp_f(1,m)=final_f(1,m)+error;
    end
%    temp_h(1,m)=final_h(1,m)+error;
%    temp_f(1,m)=final_f(1,m)+error;
   [probability_lambda,demand_lambda,~,~]=cal_pro_demand_final(delta,temp_h,final_f,velocity,sigma,v,potential_demand_true,common_true,distance_section_true); 
   [~,demand_fare,~,~]=cal_pro_demand_final(delta,final_h,temp_f,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
   partial_demand_lambda(:,m)=(demand_lambda-demand)./(temp_h(1,m)-final_h(1,m));
   partial_demand_fare(:,m)=(demand_fare-demand)./(temp_f(1,m)-final_f(1,m));
   partial_probability_lambda(:,:,m)=(probability_lambda-probability)./(temp_h(1,m)-final_h(1,m));
end
end