%基于POE模型下拉格朗日算法的机器人逆动力学求解
%目前只支持数字解
%Tor=idyna_poe_lag(robot,q,qd,qdd,gravity,f_end)
%输入robot为机器人模型，包含有运动学参数和动力学参数
%输入q为关节角度（1XN），qd为关节角加速度,qdd为对应关节角加速度
%gravity为重力加速度，f_end为末端执行器施加载荷，
%这两个在不输入时会选用默认值
%输出Torque为关节扭矩
%参考文献为：熊有伦等.机器人学建模规划与控制,7.5节:2018
%A Lie Group Formulation of Robot Dynamics

function Tor=idyna_poe_lag(robot,q,qd,qdd,gravity,f_end)

end 