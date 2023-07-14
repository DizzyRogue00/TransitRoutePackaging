%NP=10;G=5;run=10;
[transitline,transitdistance,velocity,CV,M,N,sigma,potential_demand,v,newdistance,section,total_section,~,~,~,~,max_distance,potential_demand_true,distance_section_true,common_true,section_true]=Initialization();
error=0.00001;
epsilon=0.001;
max_lambda=12;
min_lambda=1;
max_fare=10;
%max_fare=2;
min_fare=0.1;
upper=[3,3,3,3,3];
lower=[1,1,1,1,1];
[a,b,c,d,e]=ndgrid(1:3);
delta_final=[a(:),b(:),c(:),d(:),e(:)];
delta_final_row=size(delta_final,1);
delta_final_M=delta_final;
temp_eye=eye(N+1);
allowable=2;
for t=1:delta_final_row
    temp_matrix=temp_eye(delta_final(t,:),:);
    temp_sum=sum(temp_matrix,1);
    temp_result=find(temp_sum>allowable);
    if isempty(temp_result)==0
        delta_final_M(t,:)=0;
    end
end
delta_final_M(all(delta_final_M==0,2),:)=[];
delta_final_M=[(N+1).*ones(1,M);delta_final_M];
delta_final_M_row=size(delta_final_M,1);
%%
NP=10;
G=5;
run=5;
runrun=[3234,1720,1700,1630,8888];
%%

xx_M=zeros(run,M);
yy_M=zeros(run,1);
xx_M1=zeros(run,M);
yy_M1=zeros(run,1);
xx_valueM=zeros(run,2*M);
xx_valueM1=zeros(run,2*M);
xx_profitM=zeros(run,N+1);
xx_profitM1=zeros(run,N+1);
xx_LM=zeros(run,1);
xx_LM1=zeros(run,1);
xx_LLM=zeros(run,1);
xx_LLM1=zeros(run,1);
palpha=0.2;
selection=floor(NP*palpha);
for r=1:run
    rng(runrun(r));
    initial_population=randperm(delta_final_row,NP);
    initial_population_M=randperm(delta_final_M_row,NP);
    X_initialM=delta_final_M(initial_population_M,:);
    %A_initialM=X_initialM;
    %20201030
    initial_population_AM=randperm(delta_final_M_row,NP);
    A_initialM=delta_final_M(initial_population_AM,:);
    %20201030
    %%
    X_initial=delta_final(initial_population,:);
    %A_initial=X_initial;
    %20201030
    initial_population_A=randperm(delta_final_row,NP);
    A_initial=delta_final(initial_population_A,:);
    %20201030
    %%
    feasible=zeros(NP,1);
    feasible1=zeros(NP,1);
    fitness=zeros(NP,1);
    fitness1=zeros(NP,1);% service quality
    feasible_initial=zeros(NP,1);
    fitness_initial=zeros(NP,1);
    fitness_initial1=zeros(NP,1);
    %20201030
    feasible_initialA=zeros(NP,1);
    fitness_initialA=zeros(NP,1);
    fitness_initialA1=zeros(NP,1);
    %20201030
    fmin=0;
    fmin1=0;
    xbest=[3,3,3,3,3];
    xbest1=[3,3,3,3,3];
    fbest=0;
    fbest1=0;   
    %%
    %{
    20201109
    feasible_M=zeros(NP,1);
    feasible_M1=zeros(NP,1);
    fitness_M=zeros(NP,1);
    fitness_M1=zeros(NP,1);
    %}
    feasible_initialM=zeros(NP,1);
    fitness_initialM=zeros(NP,1);
    fitness_initialM1=zeros(NP,1);
    %20201030
    feasible_initialAM=zeros(NP,1);
    fitness_initialAM=zeros(NP,1);
    fitness_initialAM1=zeros(NP,1);
    %20201030
    %{
    20201109
    fmin_M=0;
    fmin_M1=0;
    xbestM=[3,3,3,3,3];
    xbestM1=[3,3,3,3,3];
    fbestM=0;
    fbestM1=0;
    %}
    %%
    tic;
    for t=1:NP
        [temp_lambda,temp_fare]=lower_model_final(X_initialM(t,:),N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
        [probability,demand,~,~]=cal_pro_demand_final(X_initialM(t,:),temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
        [profit,L]=cal_profit_socialWelfare(probability,demand,X_initialM(t,:),N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
        LL=sum(common_true.*distance_section_true.*demand.*v,'all');
        temp=find(profit(1:N)<0,1);
        fitness_initialM(t,1)=fitness_initialM(t,1)+sum(min(profit(1:N),0),2);
        fitness_initialM1(t,1)=fitness_initialM1(t,1)+sum(min(profit(1:N),0),2);
        fitness_initialM(t,1)=fitness_initialM(t,1)+min(-sum(X_initialM(t,:),2)-0.5+M*(N+1),0);
        fitness_initialM1(t,1)=fitness_initialM1(t,1)+min(-sum(X_initialM(t,:),2)-0.5+M*(N+1),0);
        num=zeros(1,N+1);
        for j=1:N+1
            num(1,j)=sum(X_initialM(t,:)==j,2);
        end
        fitness_initialM(t,1)=fitness_initialM(t,1)+sum(min(2*zeros(1,N+1)-num,0),2);
        fitness_initialM1(t,1)=fitness_initialM1(t,1)+sum(min(2*zeros(1,N+1)-num,0),2);
        if isempty(temp)
            if sum(X_initialM(t,:))<M*(N+1)
                if sum(min(2*ones(1,N+1)-num,0),2)>=0
                    feasible_initialM(t,1)=1;
                    fitness_initialM(t,1)=fitness_initialM(t,1)+L;
                    fitness_initialM1(t,1)=fitness_initialM1(t,1)+LL;
                end
            end
        end
    end
    %20201030
    for t=1:NP
        [temp_lambda,temp_fare]=lower_model_final(A_initialM(t,:),N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
        [probability,demand,~,~]=cal_pro_demand_final(A_initialM(t,:),temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
        [profit,L]=cal_profit_socialWelfare(probability,demand,A_initialM(t,:),N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
        LL=sum(common_true.*distance_section_true.*demand.*v,'all');
        temp=find(profit(1:N)<0,1);
        fitness_initialAM(t,1)=fitness_initialAM(t,1)+sum(min(profit(1:N),0),2);
        fitness_initialAM1(t,1)=fitness_initialAM1(t,1)+sum(min(profit(1:N),0),2);
        fitness_initialAM(t,1)=fitness_initialAM(t,1)+min(-sum(A_initialM(t,:),2)-0.5+M*(N+1),0);
        fitness_initialAM1(t,1)=fitness_initialAM1(t,1)+min(-sum(A_initialM(t,:),2)-0.5+M*(N+1),0);
        num=zeros(1,N+1);
        for j=1:N+1
            num(1,j)=sum(A_initialM(t,:)==j,2);
        end
        fitness_initialAM(t,1)=fitness_initialAM(t,1)+sum(min(2*zeros(1,N+1)-num,0),2);
        fitness_initialAM1(t,1)=fitness_initialAM1(t,1)+sum(min(2*zeros(1,N+1)-num,0),2);
        if isempty(temp)
            if sum(A_initialM(t,:))<M*(N+1)
                if sum(min(2*ones(1,N+1)-num,0),2)>=0
                    feasible_initialAM(t,1)=1;
                    fitness_initialAM(t,1)=fitness_initialAM(t,1)+L;
                    fitness_initialAM1(t,1)=fitness_initialAM1(t,1)+LL;
                end
            end
        end
    end
    toc;
    %20201030
    if any(feasible_initialM,'all')
        temp1=feasible_initialM~=0;
        fmin_initialM=min(fitness_initialM(temp1));
        fmin_initialM1=min(fitness_initialM1(temp1));
    else
        fmin_initialM=0;
        fmin_initialM1=0;
    end
    %20201030
    if any(feasible_initialAM,'all')
        temp1=feasible_initialAM~=0;
        fmin_initialAM=min(fitness_initialAM(temp1));
        fmin_initialAM1=min(fitness_initialAM1(temp1));
    else
        fmin_initialAM=0;
        fmin_initialAM1=0;
    end
    %20201030
    temp1=feasible_initialM==0;
    fitness_initialM(temp1)=fitness_initialM(temp1)+fmin_initialM;
    fitness_initialM1(temp1)=fitness_initialM1(temp1)+fmin_initialM1;
    %20201030
    temp1=feasible_initialAM==0;
    fitness_initialAM(temp1)=fitness_initialAM(temp1)+fmin_initialAM;
    fitness_initialAM1(temp1)=fitness_initialAM1(temp1)+fmin_initialAM1;
    %20201030
    [fbestinitialM,indM]=max(fitness_initialM);
    xbestinitialM=X_initialM(indM,:);
    [fbestinitialM1,indM1]=max(fitness_initialM1);
    xbestinitialM1=X_initialM(indM1,:);

    %%
    X_M=X_initialM;
    X_M1=X_initialM;
    A_M=A_initialM;
    A_M1=A_initialM;
    feasible_M=feasible_initialM;
    feasible_M1=feasible_initialM;
    fitness_M=fitness_initialM;
    fitness_M1=fitness_initialM1;
    fmin_M=fmin_initialM;
    fmin_M1=fmin_initialM1;
    xbestM=xbestinitialM;
    fbestM=fbestinitialM;
    xbestM1=xbestinitialM1;
    fbestM1=fbestinitialM1;
    feasible_AM=zeros(2*NP,1);
    feasible_AM1=zeros(2*NP,1);
    feasible_AM(1:NP,1)=feasible_initialAM;
    feasible_AM1(1:NP,1)=feasible_initialAM;
    fitness_AM=zeros(2*NP,1);
    fitness_AM1=zeros(2*NP,1);
    fitness_AM(1:NP,1)=fitness_initialAM;
    fitness_AM1(1:NP,1)=fitness_initialAM1;
    w_M=Dec2Binary(lower,upper,X_M);
    w_M1=Dec2Binary(lower,upper,X_M1);   
    col_num=size(w_M,2);
    w_AM=zeros(2*NP,col_num);
    w_AM1=zeros(2*NP,col_num);
    w_AM(1:NP,:)=Dec2Binary(lower,upper,A_M);
    w_AM1(1:NP,:)=Dec2Binary(lower,upper,A_M1);
    temp_M=zeros(NP,col_num);
    temp_M1=zeros(NP,col_num);
    feasible_temp=zeros(NP,1);
    feasible_temp1=zeros(NP,1);
    fitness_temp=zeros(NP,1);
    fitness_temp1=zeros(NP,1);
    %%
    scr=zeros(NP,1);
    sf=zeros(NP,1);
    F=0.5;
    CR=0.5;
    
    scr1=zeros(NP,1);
    sf1=zeros(NP,1);
    F1=0.5;
    CR1=0.5;

    c=0.2;
    c1=0.2;
    
    for gg=1:G
        tic;
        [~,order]=sort(fitness_M);
        [~,order1]=sort(fitness_M1);
        index=1;
        index1=1;
        if any(sf,'all')
            temp1=find(sf~=0);
            F=(1-c)*F+c*sum(sf(temp1).^2,1)/sum(sf(temp1),1);
        else
            if gg==1
                F=0.5;
            else
                F=(1-c)*F;
            end
        end
        if any(scr,'all')
            temp2=scr~=0;
            CR=(1-c)*CR+c*mean(scr(temp2));
        else
            if gg==1
                CR=0.5;
            else
                CR=(1-c)*CR;
            end
        end
        if any(sf1,'all')
            temp1=find(sf1~=0);
            F1=(1-c1)*F1+c1*sum(sf1(temp1).^2,1)/sum(sf1(temp1),1);
        else
            if gg==1
                F1=0.5;
            else
                F1=(1-c1)*F1;
            end
        end
        if any(scr1,'all')
            temp2=scr1~=0;
            CR1=(1-c1)*CR1+c1*mean(scr1(temp2));
        else
            if gg==1
                CR1=0.5;
            else
                CR1=(1-c1)*CR1;
            end
        end
        scr=zeros(NP,1);
        sf=zeros(NP,1);
        f_i=zeros(NP,1);
        cr_i=zeros(NP,1);
        scr1=zeros(NP,1);
        sf1=zeros(NP,1);
        f_i1=zeros(NP,1);
        cr_i1=zeros(NP,1);

        for i=1:NP
           f_i(i,1)=0.1*tan(pi*(rand()-0.5))+F;
           while 1
              if f_i(i,1)>1
                  f_i(i,1)=1;
                  break;
              end
              if f_i(i,1)>0
                  break;
              end
              f_i(i,1)=0.1*tan(pi*(rand()-0.5))+F;
           end
        end

        for i=1:NP
            cr_i(i,1)=randn()*0.1+CR;
            if cr_i(i,1)<0
                cr_i(i,1)=0;
            end
            if cr_i(i,1)>1
                cr_i(i,1)=1;
            end
        end

        for i=1:NP
           f_i1(i,1)=0.1*tan(pi*(rand()-0.5))+F1;
           while 1
              if f_i1(i,1)>1
                  f_i1(i,1)=1;
                  break;
              end
              if f_i1(i,1)>0
                  break;
              end
              f_i1(i,1)=0.1*tan(pi*(rand()-0.5))+F1;
           end
        end
        for i=1:NP
            cr_i1(i,1)=randn()*0.1+CR1;
            if cr_i1(i,1)<0
                cr_i1(i,1)=0;
            end
            if cr_i1(i,1)>1
                cr_i1(i,1)=1;
            end
        end        
        for i=1:NP
            temp_index=randperm(selection,1);
            middle_index=order(temp_index);
            middle_index1=order1(temp_index);
            wbestpM=w_M(middle_index,:);
            wbestpM1=w_M1(middle_index1,:);
            temp1=randperm(NP,1);
            while 1
                if any(w_M(temp1,:)-w_M(i,:),'all')
                    break;
                end
                temp1=randperm(NP,1);
            end
            temp2=randperm(2*NP,1);
            while 1
                if temp2<=NP
                    if any(w_M(temp2,:)-w_M(i,:),'all') && any(w_M(temp2,:)-w_M(temp1,:),'all')
                        break;
                    end
                    temp2=randperm(2*NP,1);
                else
                    temp3=temp2-NP;
                    if any(w_AM(temp3,:)-w_M(i,:),'all') && any(w_AM(temp3,:)-w_M(temp1,:),'all')
                        break;
                    end
                    temp2=randperm(2*NP,1);
                end
            end
            x1=w_M(temp1,:);
            x11=w_M1(temp1,:);
            if temp2<=NP
                x3=w_M(temp2,:);
                x33=w_M1(temp2,:);
            else
                temp2=temp2-NP;
                x3=w_AM(temp2,:);
                x33=w_AM1(temp2,:);               
            end
            v_i=w_M(i,:)+f_i(i,1)*(wbestpM-w_M(i,:))+f_i(i,1)*(x1-x3);
            v_i1=w_M1(i,:)+f_i1(i,1)*(wbestpM1-w_M1(i,:))+f_i1(i,1)*(x11-x33);
            col=col_num;
            j=randperm(col,1);
            u_i=w_M(i,:);
            for k=1:col
               if rand()<cr_i(i) || k==j
                   u_i(k)=v_i(k);
               end
            end
            u_i1=w_M(i,:);
            for k=1:col
                if rand()<cr_i1(i) || k==j
                    u_i1(k)=v_i1(k);
                end
            end
            w_i=zeros(1,col);
            w_i1=zeros(1,col);
            for k=1:col
                if u_i(k)<=0.5
                    w_i(k)=0;
                else
                    w_i(k)=1;
                end
            end
            for k=1:col
                if u_i1(k)<=0.5
                    w_i1(k)=0;
                else
                    w_i1(k)=1;
                end
            end
            w=Binary2Dec(lower,upper,w_i);
            w1=Binary2Dec(lower,upper,w_i1);
            w_num=zeros(1,N+1);
            w_num1=zeros(1,N+1);
            for jj=1:N+1
                w_num(1,jj)=sum(w==jj,2);
                w_num1(1,jj)=sum(w1==jj,2);
            end   
            feasible_temp(i,1)=0;
            if sum(min(upper-w,0),2)<0
                %randomly generate a new solution
                while 1
                    temp_i=w_i;
                    loc=find(w>upper);
                    loc_size=size(loc,2);
                    bits=ceil(log2(upper));
                    for jj=1:loc_size
                        up=sum(bits(1:loc(jj)));
                        low=up-bits(loc(jj))+1;
                        for kk=low:up
                           b=rand();
                           temp_i(kk)=b;
                           u_i(kk)=b;
                           if temp_i(kk)<=0.5
                               temp_i(kk)=0;
                           else
                               temp_i(kk)=1;
                           end
                        end
                    end
                    w_ii=Binary2Dec(lower,upper,temp_i);
                    if sum(min(upper-w_ii,0),2)>=0
                        break;
                    end
                end
                w=w_ii;
                w_num=zeros(1,N+1);
                for jj=1:N+1
                    w_num(1,jj)=sum(w==jj,2);
                end   
            end       
            [temp_lambda,temp_fare]=lower_model_final(w,N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
            [probability,demand,~,~]=cal_pro_demand_final(w,temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
            [profit,L]=cal_profit_socialWelfare(probability,demand,w,N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
            %LL=sum(common_true.*distance_section_true.*demand.*v,'all');
            temp=find(profit(1:N)<0, 1);
            if isempty(temp)
                if sum(w,2)<(N+1)*M
                    if sum(min(2*ones(1,N+1)-w_num,0),2)>=0
                        feasible_temp(i,1)=1;
                        fitness_temp(i,1)=L;
                    else
                        fitness_temp(i,1)=fmin_M+sum(min(2*ones(1,N+1)-w_num,0),2);
                    end
                else
                    fitness_temp(i,1)=fmin_M+min((N+1)*M-0.5-sum(w,2),0)+sum(min(2*ones(1,N+1)-w_num,0),2);
                end
            else
                fitness_temp(i,1)=fmin_M+sum(min(profit(1:N),0),2)+min((N+1)*M-0.5-sum(w,2),0)+sum(min(2*ones(1,N+1)-w_num,0),2);
            end  
            feasible_temp1(i,1)=0;
            if sum(min(upper-w1,0),2)<0
                %randomly generate a new solution
                while 1
                    temp_i1=w_i1;
                    loc1=find(w1>upper);
                    loc1_size=size(loc1,2);
                    bits1=ceil(log2(upper));
                    for jj=1:loc1_size
                        up=sum(bits(1:loc1(jj)));
                        low=up-bits(loc1(jj))+1;
                        for kk=low:up
                           b=rand();
                           temp_i1(kk)=b;
                           u_i1(kk)=b;
                           if temp_i1(kk)<=0.5
                               temp_i1(kk)=0;
                           else
                               temp_i1(kk)=1;
                           end
                        end
                    end
                    w_ii1=Binary2Dec(lower,upper,temp_i1);
                    if sum(min(upper-w_ii1,0),2)>=0
                        break;
                    end
                end
                w1=w_ii1;
                w_num1=zeros(1,N+1);
                for jj=1:N+1
                    w_num1(1,jj)=sum(w1==jj,2);
                end  
            end          
            [temp_lambda,temp_fare]=lower_model_final(w1,N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
            [probability,demand,~,~]=cal_pro_demand_final(w1,temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
            [profit,~]=cal_profit_socialWelfare(probability,demand,w1,N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
            LL=sum(common_true.*distance_section_true.*demand.*v,'all');
            temp=find(profit(1:N)<0);
            if isempty(temp)
                if sum(w1,2)<(N+1)*M
                    if sum(min(2*ones(1,N+1)-w_num1,0),2)>=0
                        feasible_temp1(i,1)=1;
                        fitness_temp1(i,1)=LL;
                    else
                        fitness_temp1(i,1)=fmin_M1+sum(min(2*ones(1,N+1)-w_num1,0),2);
                    end
                else
                    fitness_temp1(i,1)=fmin_M1+min((N+1)*M-0.5-sum(w1,2),0)+sum(min(2*ones(1,N+1)-w_num1,0),2);
                end
            else
                fitness_temp1(i,1)=fmin_M1+sum(min(profit(1:N),0),2)+min((N+1)*M-0.5-sum(w1,2),0)+sum(min(2*ones(1,N+1)-w_num1,0),2);
            end
            if feasible_temp(i,1)+feasible_M(i,1)~=1
                if fitness_temp(i,1)>=fitness_M(i,1)
                    temp_M(i,:)=u_i;
                    w_AM(NP+index,:)=w_M(i,:);
                    feasible_AM(NP+index,1)=feasible_M(i,:);
                    fitness_AM(NP+index,1)=fitness_M(i,:);
                    index=index+1;
                    scr(i,1)=cr_i(i,1);
                    sf(i,1)=f_i(i,1);
                else
                    temp_M(i,:)=w_M(i,:);
                    feasible_temp(i,1)=feasible_M(i,1);
                    fitness_temp(i,1)=fitness_M(i,1);
                end
            else
                if feasible_temp(i,1)==1
                    temp_M(i,:)=u_i;
                    w_AM(NP+index,:)=w_M(i,:);
                    feasible_AM(NP+index,1)=feasible_M(i,:);
                    fitness_AM(NP+index,1)=fitness_M(i,:);
                    index=index+1;
                    scr(i,1)=cr_i(i,1);
                    sf(i,1)=f_i(i,1);
                else
                    temp_M(i,:)=w_M(i,:);
                    feasible_temp(i,1)=feasible_M(i,1);
                    fitness_temp(i,1)=fitness_M(i,1);
                end
            end
            if feasible_temp1(i,1)+feasible_M1(i,1)~=1
                if fitness_temp1(i,1)>=fitness_M1(i,1)
                    temp_M1(i,:)=u_i1;
                    w_AM1(NP+index1,:)=w_M1(i,:);
                    feasible_AM1(NP+index1,1)=feasible_M1(i,:);
                    fitness_AM1(NP+index1,1)=fitness_M1(i,:);
                    index1=index1+1;
                    scr1(i,1)=cr_i1(i,1);
                    sf1(i,1)=f_i1(i,1);
                else
                    temp_M1(i,:)=w_M1(i,:);
                    feasible_temp1(i,1)=feasible_M1(i,1);
                    fitness_temp1(i,1)=fitness_M1(i,1);
                end
            else
                if feasible_temp1(i,1)==1
                    temp_M1(i,:)=u_i1;
                    w_AM1(NP+index1,:)=w_M1(i,:);
                    feasible_AM1(NP+index1,1)=feasible_M1(i,:);
                    fitness_AM1(NP+index1,1)=fitness_M1(i,:);
                    index1=index1+1;
                    scr1(i,1)=cr_i1(i,1);
                    sf1(i,1)=f_i1(i,1);
                else
                    temp_M1(i,:)=w_M1(i,:);
                    feasible_temp1(i,1)=feasible_M1(i,1);
                    fitness_temp1(i,1)=fitness_M1(i,1);
                end
            end
        end
        sizeAM=NP+index-1;
        remove_index=randperm(sizeAM,index-1);
        tempM=w_AM;
        tempM(remove_index,:)=[];
        tempM=tempM(1:NP,:);
        w_AM(1:NP,:)=tempM;
        w_AM(NP+1:2*NP,:)=zeros(NP,col_num);
        tempfeasibleM=feasible_AM;
        tempfeasibleM(remove_index,:)=[];
        tempfeasibleM=tempfeasibleM(1:NP,1);
        feasible_AM(1:NP,1)=tempfeasibleM;
        feasible_AM(NP+1:2*NP,1)=zeros(NP,1);
        tempfitnessM=fitness_AM;
        tempfitnessM(remove_index,:)=[];
        tempfitnessM=tempfitnessM(1:NP,1);
        fitness_AM(1:NP,1)=tempfitnessM;
        fitness_AM(NP+1:2*NP,1)=zeros(NP,1);
        temptemp=feasible_temp==1;
        if isempty(fitness_temp(temptemp))
            fmin_M=0;
        else
            fmin_M=min(fitness_temp(temptemp));
        end
        [fbestM,indM]=max(fitness_temp);
        xbestM=temp_M(indM,:);
        for j=1:col
            if xbestM(j)<=0.5
                xbestM(j)=0;
            else
                xbestM(j)=1;
            end
        end
        xbestM=Binary2Dec(lower,upper,xbestM);
        sizeAM1=NP+index1-1;
        remove_index1=randperm(sizeAM1,index1-1);
        tempM1=w_AM1;
        tempM1(remove_index1,:)=[];
        tempM1=tempM1(1:NP,:);
        w_AM1(1:NP,:)=tempM1;
        w_AM1(NP+1:2*NP,:)=zeros(NP,col_num);
        tempfeasibleM1=feasible_AM1;
        tempfeasibleM1(remove_index1,:)=[];
        tempfeasibleM1=tempfeasibleM1(1:NP,1);
        feasible_AM1(1:NP,1)=tempfeasibleM1;
        feasible_AM1(NP+1:2*NP,:)=zeros(NP,1);
        tempfitnessM1=fitness_AM1;
        tempfitnessM1(remove_index1,:)=[];
        tempfitnessM1=tempfitnessM1(1:NP,1);
        fitness_AM1(1:NP,1)=tempfitnessM1;
        fitness_AM1(NP+1:2*NP,1)=zeros(NP,1);
        temptemp1=feasible_temp1==1;
        if isempty(fitness_temp1(temptemp1))
            fmin_M1=0;
        else
            fmin_M1=min(fitness_temp1(temptemp1));
        end
        [fbestM1,indM1]=max(fitness_temp1);
        xbestM1=temp_M1(indM1,:);
        for j=1:col
           if xbestM1(j)<=0.5
               xbestM1(j)=0;
           else
               xbestM1(j)=1;
           end
        end
        xbestM1=Binary2Dec(lower,upper,xbestM1);
        w_M=temp_M;
        feasible_M=feasible_temp;
        fitness_M=fitness_temp;
        w_M1=temp_M1;
        feasible_M1=feasible_temp1;
        fitness_M1=fitness_temp1;
        toc
    end
    xx_M(r,:)=xbestM;
    yy_M(r,1)=fbestM;
    xx_M1(r,:)=xbestM1;
    yy_M1(r,1)=fbestM1;
    [temp_lambda,temp_fare]=lower_model_final(xbestM,N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
    [probability,demand,~,~]=cal_pro_demand_final(xbestM,temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
    [profit,L]=cal_profit_socialWelfare(probability,demand,xbestM,N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
    LL=sum(common_true.*distance_section_true.*demand.*v,'all');
    xx_valueM(r,:)=[temp_lambda,temp_fare];
    xx_profitM(r,:)=profit;
    xx_LM(r,1)=L;
    xx_LLM(r,1)=LL;
    [temp_lambda,temp_fare]=lower_model_final(xbestM1,N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
    [probability,demand,~,~]=cal_pro_demand_final(xbestM1,temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
    [profit,L]=cal_profit_socialWelfare(probability,demand,xbestM1,N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
    LL=sum(common_true.*distance_section_true.*demand.*v,'all');
    xx_valueM1(r,:)=[temp_lambda,temp_fare];
    xx_profitM1(r,:)=profit;
    xx_LM1(r,1)=L;
    xx_LLM1(r,1)=LL;
end
%%
xx=zeros(run,M);
yy=zeros(run,1);
xx1=zeros(run,M);
yy1=zeros(run,1);
xx_value=zeros(run,2*M);
xx_value1=zeros(run,2*M);
xx_profit=zeros(run,N+1);
xx_profit1=zeros(run,N+1);
xx_L=zeros(run,1);
xx_L1=zeros(run,1);
xx_LL=zeros(run,1);
xx_LL1=zeros(run,1);
for r=1:run
%%
    tic
    for t=1:NP
        [temp_lambda,temp_fare]=lower_model_final(X_initial(t,:),N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
        [probability,demand,~,~]=cal_pro_demand_final(X_initial(t,:),temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
        [profit,L]=cal_profit_socialWelfare(probability,demand,X_initial(t,:),N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
        LL=sum(common_true.*distance_section_true.*demand.*v,'all');
        temp=find(profit(1:N)<0,1);
        fitness_initial(t,1)=fitness_initial(t,1)+sum(min(profit(1:N),0),2);
        fitness_initial1(t,1)=fitness_initial1(t,1)+sum(min(profit(1:N),0),2);
        fitness_initial(t,1)=fitness_initial(t,1)+min(-sum(X_initial(t,:),2)-0.5+M*(N+1),0);
        fitness_initial1(t,1)=fitness_initial1(t,1)+min(-sum(X_initial(t,:),2)-0.5+M*(N+1),0);
        if isempty(temp)
            if sum(X_initial(t,:))<M*(N+1)
                feasible_initial(t,1)=1;
                fitness_initial(t,1)=fitness_initial(t,1)+L;
                fitness_initial1(t,1)=fitness_initial1(t,1)+LL;
            end
        end
    end
    for t=1:NP
        [temp_lambda,temp_fare]=lower_model_final(A_initial(t,:),N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
        [probability,demand,~,~]=cal_pro_demand_final(A_initial(t,:),temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
        [profit,L]=cal_profit_socialWelfare(probability,demand,A_initial(t,:),N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
        LL=sum(common_true.*distance_section_true.*demand.*v,'all');
        temp=find(profit(1:N)<0,1);
        fitness_initialA(t,1)=fitness_initialA(t,1)+sum(min(profit(1:N),0),2);
        fitness_initialA1(t,1)=fitness_initialA1(t,1)+sum(min(profit(1:N),0),2);
        fitness_initialA(t,1)=fitness_initialA(t,1)+min(-sum(A_initial(t,:),2)-0.5+M*(N+1),0);
        fitness_initialA1(t,1)=fitness_initialA1(t,1)+min(-sum(A_initial(t,:),2)-0.5+M*(N+1),0);
        if isempty(temp)
            if sum(A_initial(t,:))<M*(N+1)
                feasible_initialA(t,1)=1;
                fitness_initialA(t,1)=fitness_initialA(t,1)+L;
                fitness_initialA1(t,1)=fitness_initialA1(t,1)+LL;
            end
        end
    end
    toc
    if any(feasible_initial,'all')
        temp1=feasible_initial~=0;
        fmin_initial=min(fitness_initial(temp1));
        fmin_initial1=min(fitness_initial1(temp1));
    else
        fmin_initial=0;
        fmin_initial1=0;
    end
    if any(feasible_initialA,'all')
        temp1=feasible_initialA~=0;
        fmin_initialA=min(fitness_initialA(temp1));
        fmin_initialA1=min(fitness_initialA1(temp1));
    else
        fmin_initialA=0;
        fmin_initialA1=0;
    end
    temp1=feasible_initial==0;
    fitness_initial(temp1)=fitness_initial(temp1)+fmin_initial;
    fitness_initial1(temp1)=fitness_initial1(temp1)+fmin_initial1;
    temp1=feasible_initialA==0;
    fitness_initialA(temp1)=fitness_initialA(temp1)+fmin_initialA;
    fitness_initialA1(temp1)=fitness_initialA1(temp1)+fmin_initialA1;
    [fbestinitial,ind]=max(fitness_initial);
    xbestinitial=X_initial(ind,:);
    [fbestinitial1,ind1]=max(fitness_initial1);
    xbestinitial1=X_initial(ind1,:);

    %%
    X=X_initial;
    X1=X_initial;
    A=A_initial;
    A1=A_initial;
    feasible=feasible_initial;
    feasible1=feasible_initial;
    fitness=fitness_initial;
    fitness1=fitness_initial1;
    fmin=fmin_initial;
    fmin1=fmin_initial1;
    xbest=xbestinitial;
    fbest=fbestinitial;
    xbest1=xbestinitial1;
    fbest1=fbestinitial1;
    feasible_A=zeros(2*NP,1);
    feasible_A1=zeros(2*NP,1);
    feasible_A(1:NP,1)=feasible_initialA;
    feasible_A1(1:NP,1)=feasible_initialA;
    fitness_A=zeros(2*NP,1);
    fitness_A1=zeros(2*NP,1);
    fitness_A(1:NP,1)=fitness_initialA;
    fitness_A1(1:NP,1)=fitness_initialA1;
	w_=Dec2Binary(lower,upper,X);
	w1_=Dec2Binary(lower,upper,X1);
	col_num=size(w_,2);
	w_A=zeros(2*NP,col_num);
	w_A1=zeros(2*NP,col_num);
	w_A(1:NP,:)=Dec2Binary(lower,upper,A);
	w_A1(1:NP,:)=Dec2Binary(lower,upper,A1);
    temp_=zeros(NP,col_num);
    temp_1=zeros(NP,col_num);
    feasible__temp=zeros(NP,1);
    feasible__temp1=zeros(NP,1);
    fitness__temp=zeros(NP,1);
    fitness__temp1=zeros(NP,1);
    %%
    scr=zeros(NP,1);
    sf=zeros(NP,1);
    F=0.5;
    CR=0.5;
    scr1=zeros(NP,1);
    sf1=zeros(NP,1);
    F1=0.5;
    CR1=0.5;
    c=0.2;
    c1=0.2;
    for gg=1:G
        tic;
		[~,order]=sort(fitness);
		[~,order1]=sort(fitness1);
        index=1;
        index1=1;
        if any(sf,'all')
            temp1=find(sf~=0);
            F=(1-c)*F+c*sum(sf(temp1).^2,1)/sum(sf(temp1),1);
        else
			if gg==1
				F=0.5;
			else
				F=(1-c)*F;
			end
        end
        if any(scr,'all')
            temp2=scr~=0;
            CR=(1-c)*CR+c*mean(scr(temp2));
        else
			if gg==1
				CR=0.5;
			else
				CR=(1-c)*CR;
			end
        end
        if any(sf1,'all')
            temp1=find(sf1~=0);
            F1=(1-c1)*F1+c1*sum(sf1(temp1).^2,1)/sum(sf1(temp1),1);
        else
			if gg==1
				F1=0.5;
			else
				F1=(1-c1)*F1;
			end
        end
        if any(scr1,'all')
            temp2=scr1~=0;
            CR1=(1-c1)*CR1+c1*mean(scr1(temp2));
        else
			if gg==1
				CR1=0.5;
			else
				CR1=(1-c1)*CR1;
			end
        end
        scr=zeros(NP,1);
        sf=zeros(NP,1);
        f_i=zeros(NP,1);
        cr_i=zeros(NP,1);
        scr1=zeros(NP,1);
        sf1=zeros(NP,1);
        f_i1=zeros(NP,1);
        cr_i1=zeros(NP,1);
        for i=1:NP
           f_i(i,1)=0.1*tan(pi*(rand()-0.5))+F;
           while 1
              if f_i(i,1)>1
                  f_i(i,1)=1;
                  break;
              end
              if f_i(i,1)>0
                  break;
              end
              f_i(i,1)=0.1*tan(pi*(rand()-0.5))+F;
           end
        end
        for i=1:NP
            cr_i(i,1)=randn()*0.1+CR;
            if cr_i(i,1)<0
                cr_i(i,1)=0;
            end
            if cr_i(i,1)>1
                cr_i(i,1)=1;
            end
        end
        for i=1:NP
           f_i1(i,1)=0.1*tan(pi*(rand()-0.5))+F1;
           while 1
              if f_i1(i,1)>1
                  f_i1(i,1)=1;
                  break;
              end
              if f_i1(i,1)>0
                  break;
              end
              f_i1(i,1)=0.1*tan(pi*(rand()-0.5))+F1;
           end
        end
        for i=1:NP
            cr_i1(i,1)=randn()*0.1+CR1;
            if cr_i1(i,1)<0
                cr_i1(i,1)=0;
            end
            if cr_i1(i,1)>1
                cr_i1(i,1)=1;
            end
        end   
        for i=1:NP
			temp_index=randperm(selection,1);
			middle_index=order(temp_index);
			middle_index1=order1(temp_index);
			wbestp=w_(middle_index,:);
			wbestp1=w1_(middle_index1,:);
            temp1=randperm(NP,1);
            while 1
                if any(w_(temp1,:)-w_(i,:),'all')
                    break;
                end
                temp1=randperm(NP,1);
            end
            temp2=randperm(2*NP,1);
            while 1
				if temp2<=NP
					if any(w_(temp2,:)-w_(i,:),'all') && any(w_(temp2,:)-w_(temp1,:),'all')
						break;
					end
					temp2=randperm(2*NP,1);
				else
					temp3=temp2-NP;
					if any(w_A(temp3,:)-w_(i,:),'all') && any(w_A(temp3,:)-w_(temp1,:),'all')
						break;
					end
					temp2=randperm(2*NP,1);
				end
            end
            x1=w_(temp1,:);
            x11=w1_(temp1,:);
            if temp2<=NP
                x3=w_(temp2,:);
                x33=w1_(temp2,:);
            else
                temp2=temp2-NP;
                x3=w_A(temp2,:);
                x33=w_A1(temp2,:);
            end
			
            v_i=w_(i,:)+f_i(i,1)*(wbestp-w_(i,:))+f_i(i,1)*(x1-x3);
			v_i1=w1_(i,:)+f_i1(i,1)*(wbestp1-w1_(i,:))+f_i1(i,1)*(x11-x33);
			col=col_num;
			j=randperm(col,1);
			u_i=w_(i,:);
			for k=1:col
				if rand()<cr_i(i,1) || k==j
					u_i(k)=v_i(k);
				end
			end
			u_i1=w1_(i,:);
			for k=1:col
				if rand()<cr_i1(i,1) || k==j
					u_i1(k)=v_i1(k);
				end
			end
			w_i=zeros(1,col);
			w_i1=zeros(1,col);
			for k=1:col
				if u_i(k)<=0.5
					w_i(k)=0;
				else
					w_i(k)=1;
				end
			end
			for k=1:col
				if u_i1(k)<=0.5
					w_i1(k)=0;
				else
					w_i1(k)=1;
				end
			end
			w=Binary2Dec(lower,upper,w_i);
            w1=Binary2Dec(lower,upper,w_i1);
			w_num=zeros(1,N+1);
			w_num1=zeros(1,N+1);
			for jj=1:N+1
				w_num(1,jj)=sum(w==jj,2);
				w_num1(1,jj)=sum(w1==jj,2);
			end
            feasible__temp(i,1)=0;
            if sum(min(upper-w,0),2)<0
                %randomly generate a new solution
				while 1
					temp_i=w_i;
					loc=find(w>upper);
					loc_size=size(loc,2);
					bits=ceil(log2(upper));
					for jj=1:loc_size
						up=sum(bits(1:loc(jj)));
						low=up-bits(loc(jj))+1;
						for kk=low:up
							b=rand();
							temp_i(kk)=b;
							u_i(kk)=b;
							if temp_i(kk)<=0.5
								temp_i(kk)=0;
							else
								temp_i(kk)=1;
							end
						end
					end
					w_ii=Binary2Dec(lower,upper,temp_i);
					if sum(min(upper-w_ii,0),2)>=0
						break;
					end
				end
				w=w_ii;
				w_num=zeros(1,N+1);
				for jj=1:N+1
					w_num=sum(w==jj,2);
				end
            end
            [temp_lambda,temp_fare]=lower_model_final(w,N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
            [probability,demand,~,~]=cal_pro_demand_final(w,temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
            [profit,L]=cal_profit_socialWelfare(probability,demand,w,N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
            %LL=sum(common_true.*distance_section_true.*demand.*v,'all');
            temp=find(profit(1:N)<0, 1);
            if isempty(temp)
                if sum(w,2)<(N+1)*M
                    feasible__temp(i,1)=1;
                    fitness__temp(i,1)=L;
                else
                    fitness__temp(i,1)=fmin+min((N+1)*M-0.5-sum(w,2),0);
                end
            else
                fitness__temp(i,1)=fmin+sum(min(profit(1:N),0),2)+min((N+1)*M-0.5-sum(w,2),0);
            end
            feasible__temp1(i,1)=0;
            if sum(min(upper-w1,0),2)<0
                %randomly generate a new solution
				while 1
					temp_i1=w_i1;
					loc1=find(w1>upper);
					loc1_size=size(loc1,2);
					bits=ceil(log2(upper));
					for jj=1:loc1_size
						up=sum(bits(1:loc1(jj)));
						low=up-bits(loc1(jj))+1;
						for kk=low:up
							b=rand();
							temp_i1(kk)=b;
							u_i1(kk)=b;
							if temp_i1(kk)<=0.5
								temp_i1(kk)=0;
							else
								temp_i1(kk)=1;
							end
						end
					end
					w_ii1=Binary2Dec(lower,upper,temp_i1);
					if sum(min(upper-w_ii1,0),2)>=0
						break;
					end
				end
				w1=w_ii1;
				w_num1=zeros(1,N+1);
				for jj=1:N+1
					w_num1(1,jj)=sum(w1==jj,2);
				end
            end         
            [temp_lambda,temp_fare]=lower_model_final(w1,N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
            [probability,demand,~,~]=cal_pro_demand_final(w1,temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
            [profit,~]=cal_profit_socialWelfare(probability,demand,w1,N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
            LL=sum(common_true.*distance_section_true.*demand.*v,'all');
            temp=find(profit(1:N)<0);
            if isempty(temp)
                if sum(w1,2)<(N+1)*M
                    feasible__temp1(i,1)=1;
                    fitness__temp1(i,1)=LL;
                else
                    fitness__temp1(i,1)=fmin1+min((N+1)*M-0.5-sum(w1,2),0);
                end
            else
                fitness__temp1(i,1)=fmin1+sum(min(profit(1:N),0),2)+min((N+1)*M-0.5-sum(w1,2),0);
            end
            if feasible__temp(i,1)+feasible(i,1)~=1
                if fitness__temp(i,1)>=fitness(i,1)
                    temp_(i,:)=u_i;
                    w_A(NP+index,:)=w_(i,:);
                    feasible_A(NP+index,1)=feasible(i,:);
                    fitness_A(NP+index,1)=fitness(i,:);
                    index=index+1;
                    scr(i,1)=cr_i(i,1);
                    sf(i,1)=f_i(i,1);
                else
                    temp_(i,:)=w_(i,:);
                    feasible__temp(i,1)=feasible(i,1);
                    fitness__temp(i,1)=fitness(i,1);
                end
            else
                if feasible__temp(i,1)==1
                    temp_(i,:)=u_i;
                    w_A(NP+index,:)=w_(i,:);
                    feasible_A(NP+index,1)=feasible(i,:);
                    fitness_A(NP+index,1)=fitness(i,:);
                    index=index+1;
                    scr(i,1)=cr_i(i,1);
                    sf(i,1)=f_i(i,1);
                else
                    temp_(i,:)=w_(i,:);
                    feasible__temp(i,1)=feasible(i,1);
                    fitness__temp(i,1)=fitness(i,1);
                end
            end
            if feasible__temp1(i,1)+feasible1(i,1)~=1
                if fitness__temp1(i,1)>=fitness1(i,1)
                    temp_1(i,:)=u_i1;
                    w_A1(NP+index1,:)=w1_(i,:);
                    feasible_A1(NP+index1,1)=feasible1(i,:);
                    fitness_A1(NP+index1,1)=fitness1(i,:);
                    index1=index1+1;
                    scr1(i,1)=cr_i1(i,1);
                    sf1(i,1)=f_i1(i,1);
                else
                    temp_1(i,:)=w1_(i,:);
                    feasible__temp1(i,1)=feasible1(i,1);
                    fitness__temp1(i,1)=fitness1(i,1);
                end
            else
                if feasible__temp1(i,1)==1
                    temp_1(i,:)=u_i1;
                    w_A1(NP+index1,:)=w1_(i,:);
                    feasible_A1(NP+index1,1)=feasible1(i,:);
                    fitness_A1(NP+index1,1)=fitness1(i,:);
                    index1=index1+1;
                    scr1(i,1)=cr_i1(i,1);
                    sf1(i,1)=f_i1(i,1);
                else
                    temp_1(i,:)=w1_(i,:);
                    feasible__temp1(i,1)=feasible1(i,1);
                    fitness__temp1(i,1)=fitness1(i,1);
                end
            end
        end
        sizeA=NP+index-1;
        remove_index=randperm(sizeA,index-1);
        tempA=w_A;
        tempA(remove_index,:)=[];
        tempA=tempA(1:NP,:);
        w_A(1:NP,:)=tempA;
		w_A(NP+1:2*NP,:)=zeros(NP,col_num);
        tempfeasible=feasible_A;
        tempfeasible(remove_index,:)=[];
        tempfeasible=tempfeasible(1:NP,1);
        feasible_A(1:NP,1)=tempfeasible;
		feasible_A(NP+1:2*NP,1)=zeros(NP,1);
        tempfitness=fitness_A;
        tempfitness(remove_index,:)=[];
        tempfitness=tempfitness(1:NP,1);
        fitness_A(1:NP,1)=tempfitness;
		fitness_A(NP+1:2*NP,1)=zeros(NP,1);
        temptemp=feasible__temp==1;
        if isempty(fitness__temp(temptemp))
            fmin=0;
        else
            fmin=min(fitness__temp(temptemp));
        end
        [fbest,ind]=max(fitness__temp);
        xbest=temp_(ind,:);
		for j=1:col
			if xbest(j)<=0.5
				xbest(j)=0;
			else
				xbest(j)=1;
			end
		end
		xbest=Binary2Dec(lower,upper,xbest);
        sizeA1=NP+index1-1;
        remove_index1=randperm(sizeA1,index1-1);
        tempA1=w_A1;
        tempA1(remove_index1,:)=[];
        tempA1=tempA1(1:NP,:);
        w_A1(1:NP,:)=tempA1;
		w_A1(NP+1:2*NP,:)=zeros(NP,col_num);
        tempfeasible1=feasible_A1;
        tempfeasible1(remove_index1,:)=[];
        tempfeasible1=tempfeasible1(1:NP,1);
        feasible_A1(1:NP,1)=tempfeasible1;
		feasible_A1(NP+1:2*NP,1)=zeros(NP,1);
        tempfitness1=fitness_A1;
        tempfitness1(remove_index1,:)=[];
        tempfitness1=tempfitness1(1:NP,1);
        fitness_A1(1:NP,1)=tempfitness1;
		fitness_A1(NP+1:2*NP,1)=zeros(NP,1);
        temptemp1=feasible__temp1==1;
        if isempty(fitness__temp1(temptemp1))
            fmin1=0;
        else
            fmin1=min(fitness__temp1(temptemp1));
        end
        [fbest1,ind1]=max(fitness__temp1);
        xbest1=temp_1(ind1,:);
		for j=1:col
			if xbest1(j)<=0.5
				xbest1(j)=0;
			else
				xbest1(j)=1;
			end
		end
		xbest1=Binary2Dec(lower,upper,xbest1);
        w_=temp_;
        feasible=feasible__temp;
        fitness=fitness__temp;
        w1_=temp_1;
        feasible1=feasible__temp1;
        fitness1=fitness__temp1;
        toc
    end
    xx(r,:)=xbest;
    yy(r,1)=fbest;
    xx1(r,:)=xbest1;
    yy1(r,1)=fbest1;
    [temp_lambda,temp_fare]=lower_model_final(xbest,N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
    [probability,demand,~,~]=cal_pro_demand_final(xbest,temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
    [profit,L]=cal_profit_socialWelfare(probability,demand,xbest,N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
    LL=sum(common_true.*distance_section_true.*demand.*v,'all');
    xx_value(r,:)=[temp_lambda,temp_fare];
    xx_profit(r,:)=profit;
    xx_L(r,1)=L;
    xx_LL(r,1)=LL;
    [temp_lambda,temp_fare]=lower_model_final(xbest1,N,M,velocity,sigma,v,potential_demand_true,common_true,distance_section_true,max_distance,CV);
    [probability,demand,~,~]=cal_pro_demand_final(xbest1,temp_lambda,temp_fare,velocity,sigma,v,potential_demand_true,common_true,distance_section_true);
    [profit,L]=cal_profit_socialWelfare(probability,demand,xbest1,N,common_true,distance_section_true,temp_fare,temp_lambda,max_distance,velocity,CV,v);
    LL=sum(common_true.*distance_section_true.*demand.*v,'all');
    xx_value1(r,:)=[temp_lambda,temp_fare];
    xx_profit1(r,:)=profit;
    xx_L1(r,1)=L;
    xx_LL1(r,1)=LL;
end

%save('final_evolution_maxfare_10_20201109.mat')
save('final_evolution_maxfare_2_20201110.mat')