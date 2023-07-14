function [profit,L]=cal_profit_socialWelfare(probability,demand,delta,N,common_true,distance_section_true,final_f,final_h,max_distance,velocity,CV,v)
profit=zeros(1,N+1);
for n=1:N+1
    line_set=find(delta==n);
    if isempty(line_set)==0
        for m=1:size(line_set,2)
            index_line=find(common_true(:,line_set(m))>0);
            profit(1,n)=profit(1,n)+sum(probability(index_line,line_set(m)).*demand(index_line).*distance_section_true(index_line).*final_f(line_set(m)),1)-max_distance(line_set(m))/velocity(1,n)*final_h(line_set(m))*CV(1,n);
            %profit(1,n)=profit(1,n)+sum(probability(index_line,line_set(m)).*demand(index_line).*distance_section_true(index_line).*final_f(line_set(m)),1);
        end
%         if n==N+1
%             L=sum(common_true.*distance_section_true.*demand.*v,'all')+sum(profit,2);
%         end
    end
end
L=sum(common_true.*distance_section_true.*demand.*v,'all')+sum(profit,2);
end