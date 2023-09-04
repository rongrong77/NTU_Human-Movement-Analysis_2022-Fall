import numpy as np
import matplotlib.pyplot as plt
from Function import*

BW = 58  # 體重
vicon_data, vicon_data_name, frame_rate, analog_rate = readc3d('walking.c3d')  # 匯入 .c3d 檔資料
vicon_data_1, vicon_data_name_1, frame_rate_1, analog_rate_1 = readc3d('subcali.c3d')
start_frame, end_frame = gait_cycle(BW, vicon_data, frame_rate, analog_rate)  # 計算步態週期的起始幀和結束幀


#pelvis
RASI = np.array(vicon_data.get('RASI'))[int(start_frame):int(end_frame)+1,0:3]
LASI = np.array(vicon_data.get('LASI'))[int(start_frame):int(end_frame)+1,0:3]
RPSI = np.array(vicon_data.get('RPSI'))[int(start_frame):int(end_frame)+1,0:3]

rPcoord=np.empty([3,169,3])
rPcoord=np.array([RASI,LASI,RPSI])
lPcoord=np.empty([3,169,3])
lPcoord=np.array([RASI,LASI,RPSI])

rRg2p,rVg2p,rPcoord_local = CoordPelvis(rPcoord, 'r')
lRg2p,lVg2p,lPcoord_local = CoordPelvis(lPcoord, 'l')

#thigh-r
RTRO = np.array(vicon_data.get('RTRO'))[int(start_frame):int(end_frame)+1,0:3]
RLFC = np.array(vicon_data.get('RLFC'))[int(start_frame):int(end_frame)+1,0:3]
RMFC = np.array(vicon_data.get('RMFC'))[int(start_frame):int(end_frame)+1,0:3]

rTcoord=np.empty([3,169,3])
rTcoord=np.array([RTRO,RLFC,RMFC])

rRg2t,rVg2t,rTcoord_local = CoordThigh(rTcoord, 'r')

#thigh-l
LTRO = np.array(vicon_data.get('LTRO'))[int(start_frame):int(end_frame)+1,0:3]
LLFC = np.array(vicon_data.get('LLFC'))[int(start_frame):int(end_frame)+1,0:3]
LMFC = np.array(vicon_data.get('LMFC'))[int(start_frame):int(end_frame)+1,0:3]

lTcoord=np.empty([3,169,3])
lTcoord=np.array([LTRO,LLFC,LMFC])

lRg2t,lVg2t,lTcoord_local = CoordThigh(lTcoord, 'l')


#shank-r
RTT = np.array(vicon_data.get('RTT'))[int(start_frame):int(end_frame)+1,0:3]
RSHA = np.array(vicon_data.get('RSHA'))[int(start_frame):int(end_frame)+1,0:3]
RLMA = np.array(vicon_data.get('RLMA'))[int(start_frame):int(end_frame)+1,0:3]
RMMA= np.array(vicon_data.get('RMMA'))[int(start_frame):int(end_frame)+1,0:3]

rScoord=np.empty([4,169,3])
rScoord = np.array([RTT,RSHA,RLMA,RMMA])

rRg2s,rVg2s,rScoord_local = CoordShank(rScoord, 'r')

#shank-l
LTT = np.array(vicon_data.get('LTT'))[int(start_frame):int(end_frame)+1,0:3]
LSHA = np.array(vicon_data.get('LSHA'))[int(start_frame):int(end_frame)+1,0:3]
LLMA = np.array(vicon_data.get('LLMA'))[int(start_frame):int(end_frame)+1,0:3]
LMMA= np.array(vicon_data.get('LMMA'))[int(start_frame):int(end_frame)+1,0:3]

lScoord=np.empty([4,169,3])
lScoord = np.array([LTT,LSHA,LLMA,LMMA])

lRg2s,lVg2s,lScoord_local = CoordShank(lScoord, 'l')


#foot-r
RHEE = np.array(vicon_data.get('RHEE'))[int(start_frame):int(end_frame)+1,0:3]
RFOO = np.array(vicon_data.get('RFOO'))[int(start_frame):int(end_frame)+1,0:3]
RTOE = np.array(vicon_data.get('RTOE'))[int(start_frame):int(end_frame)+1,0:3]

rFcoord = np.array([RHEE,RFOO,RTOE])

rRg2f,rVg2f,rFcoord_local = CoordFoot(rFcoord, 'r')

LHEE = np.array(vicon_data.get('LHEE'))[int(start_frame):int(end_frame)+1,0:3]
LFOO = np.array(vicon_data.get('LFOO'))[int(start_frame):int(end_frame)+1,0:3]
LTOE = np.array(vicon_data.get('LTOE'))[int(start_frame):int(end_frame)+1,0:3]

lFcoord = np.array([LHEE,LFOO,LTOE])
lRg2f,lVg2f,lFcoord_local = CoordFoot(lFcoord, 'l')

#1

theta_lRg2p = Rot2Ang(lRg2p, 'zxy')
theta_rRg2p = Rot2Ang(rRg2p, 'zxy')

xi = theta_lRg2p[:,2]
n=169
yi = range(0,n)

p = np.polyfit(xi, yi, 5) 
dp_1 = np.polyder(p,1)
dp_2 = np.polyder(p,2)
yi_1 = np.polyval(p, xi)

smprate = 120
ltAV= Ang2LocalAngular(Rot2Ang(lRg2t,'zxy'), 'zxy', smprate )




'''
下方程式碼將把上方計算之變數存為 .npy 檔，以利批改。
'''


