function [profit,L]=cal_profit(probability,demand,delta,M,N,common_true,distance_section_true,final_f,final_h,max_distance,velocity,CV,dual_multi,min_h,max_h,min_f,max_f,mu,v)
profit=zeros(1,N+1);
L=zeros(1,N+1);
for n=1:N+1
    line_set=find(delta==n);
    if isempty(line_set)==0
        for m=1:size(line_set,2)
            index_line=find(common_true(:,line_set(m))>0);
            profit(1,n)=profit(1,n)+sum(probability(index_line,line_set(m)).*demand(index_line).*distance_section_true(index_line).*final_f(line_set(m)),1)-max_distance(line_set(m))/velocity(1,n)*final_h(line_set(m))*CV(1,n);
            %profit(1,n)=profit(1,n)+sum(probability(index_line,line_set(m)).*demand(index_line).*distance_section_true(index_line).*final_f(line_set(m)),1);
        end
        line_num=size(line_set,2);
        if n~=N+1
            dual=zeros(1,4*line_num);
            con=zeros(1,4*line_num);
        else
            dual=zeros(1,4*line_num+1);
            con=zeros(1,4*line_num+1);
        end
        for item=1:2:2*line_num-1
            dual(item)=dual_multi(2*line_set((item+1)/2)-1);
            con(item)=final_h(1,line_set((item+1)/2))-min_h;
            dual(item+1)=dual_multi(2*line_set((item+1)/2));
            con(item+1)=max_h-final_h(1,line_set((item+1)/2));
            dual(2*line_num+item)=dual_multi(2*M+2*line_set((item+1)/2)-1);
            con(2*line_num+item)=final_f(1,line_set((item+1)/2))-min_f;
            dual(2*line_num+item+1)=dual_multi(2*M+2*line_set((item+1)/2));
            con(2*line_num+item+1)=max_f-final_f(1,line_set((item+1)/2));
        end
        if n==N+1
            dual(4*line_num+1)=dual_multi(4*M+1);
            con(4*line_num+1)=profit(1,n);
        end
        if n~=N+1
            L(1,n)=-profit(1,n)+1/2/mu(n)*sum(max(0,dual-mu(1,n)*con).^2-dual.^2,2);
        else
            L(1,n)=-sum(common_true.*distance_section_true.*demand.*v,'all')-sum(profit,2)+1/2/mu(n)*sum(max(0,dual-mu(n)*con).^2-dual.^2,2);
        end
    end
end
end