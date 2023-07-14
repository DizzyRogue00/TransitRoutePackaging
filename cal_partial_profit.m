%function [sum_partial_profit_lambda,sum_partial_profit_fare]=cal_partial_profit(probability,demand,partial_demand_lambda,partial_demand_fare,partial_probability_lambda,delta,M,N,final_h,final_f,velocity,CV,v,potential_demand_true,common_true,distance_section_true,max_distance,dual_multi,mu,max_h,min_h,max_f,min_f,sigma)
function [sum_partial_profit_lambda,sum_partial_profit_fare]=cal_partial_profit(probability,demand,partial_demand_lambda,partial_demand_fare,partial_probability_lambda,delta,M,N,final_h,final_f,velocity,CV,v,common_true,distance_section_true,max_distance,dual_multi,mu,max_h,min_h,max_f,min_f)
sum_partial_profit_lambda=zeros(1,M);
sum_partial_profit_fare=zeros(1,M);
%max_lambda=max_h*ones(1,M);
%max_fare=max_f*ones(1,M);
for n=1:N
    line_set=find(delta==n);
    if isempty(line_set)==0
        for m=1:size(line_set,2)
            for m_=1:size(line_set,2)
                if m==m_
                    index_section=find(common_true(:,line_set(m))>0);
                    sum_partial_profit_lambda(1,line_set(m))=sum_partial_profit_lambda(1,line_set(m))-(sum(partial_probability_lambda(index_section,line_set(m),line_set(m)).*demand(index_section).*distance_section_true(index_section).*final_f(1,line_set(m))+...
                        probability(index_section,line_set(m)).*partial_demand_lambda(index_section,line_set(m)).*distance_section_true(index_section).*final_f(1,line_set(m)),1)-...
                        max_distance(1,line_set(m))./velocity(1,n).*CV(1,n));
                    sum_partial_profit_fare(1,line_set(m))=sum_partial_profit_fare(1,line_set(m))-sum(probability(index_section,line_set(m)).*demand(index_section).*distance_section_true(index_section)+...
                        probability(index_section,line_set(m)).*partial_demand_fare(index_section,line_set(m)).*distance_section_true(index_section).*final_f(1,line_set(m)),1);
                else
                    index_section=find(common_true(:,line_set(m_))>0);
                    sum_partial_profit_lambda(1,line_set(m))=sum_partial_profit_lambda(1,line_set(m))-sum(partial_probability_lambda(index_section,line_set(m_),line_set(m)).*demand(index_section).*distance_section_true(index_section).*final_f(1,line_set(m_))+...
                        probability(index_section,line_set(m_)).*partial_demand_lambda(index_section,line_set(m)).*distance_section_true(index_section).*final_f(1,line_set(m_)),1);
                    sum_partial_profit_fare(1,line_set(m))=sum_partial_profit_fare(1,line_set(m))-sum(probability(index_section,line_set(m_)).*partial_demand_fare(index_section,line_set(m)).*distance_section_true(index_section).*final_f(1,line_set(m_)),1);
                end
            end
            sum_partial_profit_lambda(1,line_set(m))=sum_partial_profit_lambda(1,line_set(m))-max(0,dual_multi(2*line_set(m)-1:2*line_set(m))-mu(n).*[final_h(1,line_set(m))-min_h,max_h-final_h(1,line_set(m))])*[1;-1];
            %sum_partial_profit_lambda(1,line_set(m))=-sum_partial_profit_lambda(1,line_set(m));
            sum_partial_profit_fare(1,line_set(m))=sum_partial_profit_fare(1,line_set(m))-max(0,dual_multi(2*M+2*line_set(m)-1:2*M+2*line_set(m))-mu(n).*[final_f(1,line_set(m))-min_f,max_f-final_f(1,line_set(m))])*[1;-1];
        end
    end
end
line_set=find(delta==(N+1));
% if isempty(line_set)==0
%     [~,L]=cal_profit(probability,demand,delta,M,N,common_true,distance_section_true,final_f,final_h,max_distance,velocity,CV,dual_multi,min_h,max_h,min_f,max_f,mu,v);
%     %item2=sum(profit(1:N));
%     %item2=sum(profit);
%     %construct function
%     %L=-sum(common_true.*distance_section_true.*demand.*v,'all')-item2-dual_multi*profit(N+1)+1/2*mu*(profit(N+1)^2);
%     %L=-sum(common_true.*distance_section_true.*demand.*v,'all')-item2+1/2/mu*(max(0,dual_multi-mu*profit(N+1))^2-dual_multi^2);
%     error=0.0001;
%     epsilon=0.0001;
%     for m=1:size(line_set,2)
%         temp_h=final_h;
%         temp_f=final_f;
%         if abs(final_h(1,line_set(m))-max_lambda(1,line_set(m)))<=epsilon
%             temp_h(1,line_set(m))=final_h(1,line_set(m))-error;
%         else
%             temp_h(1,line_set(m))=final_h(1,line_set(m))+error;
%         end
%         if abs(final_f(1,line_set(m))-max_fare(1,line_set(m)))<=epsilon
%             temp_f(1,line_set(m))=final_f(1,line_set(m))-error;
%         else
%             temp_f(1,line_set(m))=final_f(1,line_set(m))+error;
%         end
%         [probability_h,demand_h,~,~]=cal_pro_demand_final(delta,temp_h,final_f,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
%         [~,L_h]=cal_profit(probability_h,demand_h,delta,M,N,common_true,distance_section_true,final_f,temp_h,max_distance,velocity,CV,dual_multi,min_h,max_h,min_f,max_f,mu,v);
%         %item2_h=sum(profit_h(1:N));
%         %item2_h=sum(profit_h);
%         %L_h=-sum(common_true.*distance_section_true.*demand_h.*v,'all')-item2_h-dual_multi*profit_h(N+1)+1/2*mu*(profit_h(N+1)^2);
%         %L_h=-sum(common_true.*distance_section_true.*demand_h.*v,'all')-item2_h+1/2/mu*(max(0,dual_multi-mu*profit_h(N+1))^2-dual_multi^2);
%         sum_partial_profit_lambda(1,line_set(m))=(L_h(N+1)-L(N+1))./(temp_h(line_set(m))-final_h(line_set(m)));
%         [probability_f,demand_f,~,~]=cal_pro_demand_final(delta,final_h,temp_f,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
%         [~,L_f]=cal_profit(probability_f,demand_f,delta,M,N,common_true,distance_section_true,temp_f,final_h,max_distance,velocity,CV,dual_multi,min_h,max_h,min_f,max_f,mu,v);
%         %item2_f=sum(profit_f(1:N));
%         %item2_f=sum(profit_f);
%         %L_f=-sum(common_true.*distance_section_true.*demand_f.*v,'all')-item2_f-dual_multi*profit_f(N+1)+1/2*mu*(profit_f(N+1)^2);
%         %L_f=-sum(common_true.*distance_section_true.*demand_f.*v,'all')-item2_f+1/2/mu*(max(0,dual_multi-mu*profit_f(N+1))^2-dual_multi^2);
%         sum_partial_profit_fare(1,line_set(m))=(L_f(N+1)-L(N+1))./(temp_f(line_set(m))-final_f(line_set(m)));
%     end
% end
%% %2020/10/02
if isempty(line_set)==0
    [profit,~]=cal_profit(probability,demand,delta,M,N,common_true,distance_section_true,final_f,final_h,max_distance,velocity,CV,dual_multi,min_h,max_h,min_f,max_f,mu,v);
    for m=1:size(line_set,2)
        a1=-sum(common_true.*partial_demand_lambda(:,line_set(m)).*distance_section_true.*v,'all');
        a2=-sum(common_true.*partial_demand_fare(:,line_set(m)).*distance_section_true.*v,'all');
        for m_=1:size(line_set,2)
            if m==m_
                index_section=find(common_true(:,line_set(m))>0);
                sum_partial_profit_lambda(1,line_set(m))=sum_partial_profit_lambda(1,line_set(m))-(sum(partial_probability_lambda(index_section,line_set(m),line_set(m)).*demand(index_section).*distance_section_true(index_section).*final_f(1,line_set(m))+...
                    probability(index_section,line_set(m)).*partial_demand_lambda(index_section,line_set(m)).*distance_section_true(index_section).*final_f(1,line_set(m)),1)-...
                    max_distance(1,line_set(m))./velocity(1,N+1).*CV(1,N+1));
                sum_partial_profit_fare(1,line_set(m))=sum_partial_profit_fare(1,line_set(m))-sum(probability(index_section,line_set(m)).*demand(index_section).*distance_section_true(index_section)+...
                    probability(index_section,line_set(m)).*partial_demand_fare(index_section,line_set(m)).*distance_section_true(index_section).*final_f(1,line_set(m)),1);
            else
                index_section=find(common_true(:,line_set(m_))>0);
                sum_partial_profit_lambda(1,line_set(m))=sum_partial_profit_lambda(1,line_set(m))-sum(partial_probability_lambda(index_section,line_set(m_),line_set(m)).*demand(index_section).*distance_section_true(index_section).*final_f(1,line_set(m_))+...
                    probability(index_section,line_set(m_)).*partial_demand_lambda(index_section,line_set(m)).*distance_section_true(index_section).*final_f(1,line_set(m_)),1);
                sum_partial_profit_fare(1,line_set(m))=sum_partial_profit_fare(1,line_set(m))-sum(probability(index_section,line_set(m_)).*partial_demand_fare(index_section,line_set(m)).*distance_section_true(index_section).*final_f(1,line_set(m_)),1);
            end
        end   
        temp_1=sum_partial_profit_lambda(1,line_set(m));
        temp_2=sum_partial_profit_fare(1,line_set(m));
        sum_partial_profit_lambda(1,line_set(m))=sum_partial_profit_lambda(1,line_set(m))-max(0,dual_multi(4*M+1)-mu(N+1)*profit(N+1))*(-temp_1);
        sum_partial_profit_fare(1,line_set(m))=sum_partial_profit_fare(1,line_set(m))-max(0,dual_multi(4*M+1)-mu(N+1)*profit(N+1))*(-temp_2);
        sum_partial_profit_lambda(1,line_set(m))=sum_partial_profit_lambda(1,line_set(m))+a1;
        sum_partial_profit_fare(1,line_set(m))=sum_partial_profit_fare(1,line_set(m))+a2;
        line_set_no=find(delta~=(N+1));
        for m_=1:size(line_set_no,2)
            index_section=find(common_true(:,line_set_no(m_))>0);
            sum_partial_profit_lambda(1,line_set(m))=sum_partial_profit_lambda(1,line_set(m))-sum(partial_probability_lambda(index_section,line_set_no(m_),line_set(m)).*demand(index_section).*distance_section_true(index_section).*final_f(1,line_set_no(m_))+...
                probability(index_section,line_set_no(m_)).*partial_demand_lambda(index_section,line_set(m)).*distance_section_true(index_section).*final_f(1,line_set_no(m_)),1);
            sum_partial_profit_fare(1,line_set(m))=sum_partial_profit_fare(1,line_set(m))-sum(probability(index_section,line_set_no(m_)).*partial_demand_fare(index_section,line_set(m)).*distance_section_true(index_section).*final_f(1,line_set_no(m_)),1);
        end
        sum_partial_profit_lambda(1,line_set(m))=sum_partial_profit_lambda(1,line_set(m))-max(0,dual_multi(2*line_set(m)-1:2*line_set(m))-mu(N+1).*[final_h(1,line_set(m))-min_h,max_h-final_h(1,line_set(m))])*[1;-1];
        sum_partial_profit_fare(1,line_set(m))=sum_partial_profit_fare(1,line_set(m))-max(0,dual_multi(2*M+2*line_set(m)-1:2*M+2*line_set(m))-mu(N+1).*[final_f(1,line_set(m))-min_f,max_f-final_f(1,line_set(m))])*[1;-1];
    end
end
end