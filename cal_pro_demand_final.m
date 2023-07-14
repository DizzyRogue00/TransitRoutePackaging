function [probability,demand,EW,ET]=cal_pro_demand_final(delta,final_h,final_f,velocity,sigma,v,potential_demand_true,common_true,distance_section_true)
velocity_=velocity(delta);
probability=final_h.*common_true./sum(common_true.*final_h,2);
EW=1./sum(common_true.*final_h,2);
ES=sum((distance_section_true./velocity_).*common_true.*final_h,2)./sum(common_true.*final_h,2);
ET=EW+ES;
demand=potential_demand_true.*exp(-(sum(final_f.*distance_section_true.*probability,2)+ET.*sigma)./v./distance_section_true);
end