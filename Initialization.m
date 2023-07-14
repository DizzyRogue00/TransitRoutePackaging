function [transitline,transitdistance,velocity,CV,M,N,sigma,potential_demand,v,newdistance,section,total_section,distance_section,common,potential_demand_plus,potential_demand_matrix,max_distance,potential_demand_true,distance_section_true,common_true,section_true]=Initialization()
    %line1=[1,2,5,7,9,7,5,2,1];
    %distance1=[0,6,8,10,14,18,20,22,28];
    line1=[1,2,5,7,8,7,5,2,1];
    distance1=[0,6,8,10,12,14,16,18,24];
    line2=[3,2,5,7,8,9,8,7,5,2,3];
    distance2=[0,4,6,8,10,12,14,16,18,20,24];
    line3=[4,5,8,9,8,5,4];
    distance3=[0,4,10,12,14,20,24];
    line4=[6,7,8,7,6];
    distance4=[0,8,10,12,20];
    line5=[9,10,9];
    distance5=[0,8,16];
    velocity=[30,15,20];
    M=5;
    N=2;
    CV=[240,105,150];
    sigma=60;
    potential_demand1=[
        0,150,0,0,300,0,150,100,0,0;
        0,0,150,0,250,0,200,300,200,0;
        0,0,0,0,150,0,200,250,100,0;
        0,0,0,0,200,0,0,250,150,0;
        0,0,0,0,0,0,250,450,150,0;
        0,0,0,0,0,0,200,250,0,0;
        0,0,0,0,0,0,0,350,100,0;
        0,0,0,0,0,0,0,0,150,0;
        0,0,0,0,0,0,0,0,0,50;
        0,0,0,0,0,0,0,0,0,0
    ];
    potential_demand=potential_demand1'+potential_demand1;
    v=3;
    %v=1;
    transitline=cell(1,M);
    transitdistance=cell(1,M);
    for i=1:M
        transitline{i}=eval(['line',num2str(i)]);
        transitdistance{i}=eval(['distance',num2str(i)]);
    end
    newdistance=cell(1,M);
    for i=1:M
        transitcol=size(transitline{i},2);
        mid=(transitcol+1)/2;
        newdistance{i}=zeros(mid,mid);
        for j=1:mid-1
            for k=j+1:mid
                newdistance{i}(j,k)=transitdistance{i}(1,k)-transitdistance{i}(1,j);
            end
        end
        for j=mid:transitcol
            for k=j+1:transitcol
                jj=2*mid-j;
                kk=2*mid-k;
                newdistance{i}(jj,kk)=transitdistance{i}(1,k)-transitdistance{i}(1,j);
            end
        end
    end
    section=zeros(1,M);
    for i=1:M
        transitcol=size(transitline{i},2);
        mid=(transitcol+1)/2-1;
        section(1,i)=(1+mid)*mid;
    end
    total_section=sum(section);
    distance_section=zeros(1,total_section);
    temp_sum=0;
    for i=1:M
        transitcol=size(transitline{i},2);
        mid=(transitcol+1)/2;
        temp=reshape(newdistance{i},[1,mid*mid]);
        temp(temp==0)=[];
        if i==1
            distance_section(1,1:section(1,i))=temp;
        else
            distance_section(1,temp_sum+1:temp_sum+section(1,i))=temp;
        end
        temp_sum=temp_sum+section(1,i);
    end
    common=zeros(total_section,M);
    newcommon=cell(1,M);
    for i=1:M
        newcommon{i}=cell(1,M);
        for j=1:M
            transitcol=size(transitline{j},2);
            mid=(transitcol+1)/2;
            newcommon{i}{j}=zeros(mid,mid);  
            for k=1:mid-1
                for p=k+1:mid
                    if ismember(transitline{j}(1,k),transitline{i}) && ismember(transitline{j}(1,p),transitline{i})
                        newcommon{i}{j}(k,p)=1;
                    end
                end
            end
            for k=mid:transitcol-1
                for p=k+1:transitcol
                    if ismember(transitline{j}(1,k),transitline{i}) && ismember(transitline{j}(1,p),transitline{i})
                       kk=2*mid-k;
                       pp=2*mid-p;
                       newcommon{i}{j}(kk,pp)=1;
                    end
                end
            end
        end
    end
    for i=1:M
        temp=[];
        for j=1:M
            transitcol=size(transitline{j},2);
            mid=(1+transitcol)/2;
            temp_temp=1:(mid+1):mid*mid;
            temp1=reshape(newcommon{i}{j}.',1,[]);
            temp1(temp_temp)=[];
            temp=[temp,temp1];
        end
        common(:,i)=temp.';
    end
    potential_demand_plus=zeros(total_section,1);
    potential_demand_matrix=cell(1,M);
    t_sum=0;
    for i=1:M
        transitcol=size(transitline{i},2);
        mid=(transitcol+1)/2;
        potential_demand_matrix{i}=zeros(mid,mid);
        for j=1:mid-1
            for k=j+1:mid
                potential_demand_matrix{i}(j,k)=potential_demand(transitline{i}(1,j),transitline{i}(1,k));
            end
        end
        for j=mid:transitcol-1
            for k=j+1:transitcol
                jj=2*mid-j;
                kk=2*mid-k;
                potential_demand_matrix{i}(jj,kk)=potential_demand(transitline{i}(1,j),transitline{i}(1,k));
            end
        end
        index=1:(mid+1):mid*mid;
        index_1=reshape(potential_demand_matrix{i}.',1,[]);
        index_1(index)=[];
        potential_demand_plus(t_sum+1:t_sum+size(index_1,2),1)=index_1.';
        t_sum=t_sum+size(index_1,2);
    end
    max_distance=zeros(1,M);
    for i=1:M
        transitcol=size(transitline{i},2);
        max_distance(1,i)=transitdistance{i}(1,transitcol);
    end
    %total_section,distance_section,common,potential_demand_plus
    potential_demand_true=reshape(potential_demand,[],1);
    [row,col]=size(potential_demand);
    distance_section_true=zeros(row*col,1);
    common_true=zeros(row*col,M);
    for i=1:row*col
        %row_num=mod(i-1,row)+1;
        %col_num=floor((i-1)/row)+1;
        row_num=floor((i-1)/col)+1;
        col_num=mod(i-1,col)+1;
        for j=1:M
            if ismember(row_num,transitline{j}) && ismember(col_num,transitline{j})
                true_1=find(transitline{j}==row_num,1);
                true_2=find(transitline{j}==col_num,1);
                distance_section_true(i,1)=abs(transitdistance{j}(1,true_2)-transitdistance{j}(1,true_1));
                common_true(i,j)=1;
            end
        end
    end
    index_true=find(potential_demand_true==0);
    potential_demand_true(potential_demand_true==0)=[];
    distance_section_true(index_true)=[];
    common_true(index_true,:)=[];
    section_true=size(potential_demand_true,1);
end