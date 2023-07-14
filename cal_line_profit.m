function [line_profit]=cal_line_profit(probability,demand,delta,distance_section_true,final_f,final_h,max_distance,velocity,CV)
line_profit=sum(probability.*demand.*distance_section_true.*final_f,1)-max_distance./velocity(delta).*final_h.*CV(delta);
end