

system_setting.AOI = [750 1100; 1250 1600];                                 % Area of interest
sel_veh_id  = 20;                                                           % Selected vehicle ID
viz_frame_rate = 0.1;                                                       % Frame rate for visualization
dt = 0.05;

simulation_scenario = 3;

if simulation_scenario == 1
    
    system_setting.sigma_GPS = 5.0;
    system_setting.sigma_v = 0.5;
    system_setting.sigma_IMU = (5/180*pi);
    system_setting.sigma_UWB = 0.3;
    system_setting.com_radius = 200;
    system_setting.dt = dt;
    gps_agent = 1:250;
    infra_nodes = [];

elseif simulation_scenario == 2
    
    system_setting.sigma_GPS = 0.1;
    system_setting.sigma_v = 0.5;
    system_setting.sigma_IMU = (5/180*pi);
    system_setting.sigma_UWB = 0.3;
    system_setting.com_radius = 100;
    system_setting.dt = dt;

    %gps_agent = [13 10];
    gps_agent = [];

    %infra_nodes = [1002 1000 1300; 1003 1050 1400; 1004 1050 1500; 1005 1000 1550; 1006 1450 1000; 1007 900 1570; 1009 800 1540; 1010 900 1280; 1011 800 1280; 1012 1000 1450; 1013 1100 1280];
    infra_nodes = [1002 1000 1300; 1003 1050 1400; 1005 1000 1550; 1006 1450 1000; 1009 800 1540; 1011 800 1280; 1012 1000 1450; 1013 1100 1280];
    
elseif simulation_scenario == 3
    
    system_setting.sigma_GPS = 0.1;
    system_setting.sigma_v = 5;
    system_setting.sigma_IMU = (15/180*pi);
    system_setting.sigma_UWB = 0.3;
    system_setting.com_radius = 200;
    system_setting.dt = dt;

    gps_agent = [];
    infra_nodes = [1001 800 1540; 1002 800 1280; 1003 1000 1300; 1005 1050 1400];
end
