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
    CV=[240,105];
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
    transitline=cell(1,M);
    transitdistance=cell(1,M);
    for i=1:M
        transitline{i}=eval(['line',num2str(i)]);
        transitdistance{i}=eval(['distance',num2str(i)]);
    end
    newdistance=cell(1,M);
    

end