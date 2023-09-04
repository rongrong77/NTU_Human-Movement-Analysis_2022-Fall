import numpy as np
from Function import*

BW = 58  # 體重
vicon_data, vicon_data_name, frame_rate, analog_rate = readc3d('walking.c3d')  # 匯入 .c3d 檔資料
start_frame, end_frame = gait_cycle(BW, vicon_data, frame_rate, analog_rate)  # 計算步態週期的起始幀和結束幀


'''
從這裡開始為各位作業一需撰寫的主程式部分
以下將簡要講解上面計算出來的變數（以課堂講解為主）：

1. vicon_data：為 .c3d 檔中的 markers 與力板資料，由 dictionary 的形式儲存。力板的
資料結構為一維的 list，大小為 [nframes]，取樣頻率為 analog_rate。markers 的資料結構
為二維的 list，大小為 [nframes×5]，取樣頻率為 frame_rate。
2. vicon_data_name：為 .c3d 檔中的變數名稱，包含 vicon_data 中的所有 key，可用於提
取 vicon_data 的 markers 與力板資料。
3. frame_rate：為 .c3d 檔中 markers 的取樣頻率。
4. analog_rate：為 .c3d 檔中力板的取樣頻率。
5. start_frame：為步態週期的起始幀，可用來取得一步態週期開始時各 marker 的位置向量。
6. end_frame：為步態週期的結束幀，可用來取得一步態週期結束時各 marker 的位置向量。


請於下方運用各位在 function.py 撰寫的函式，修改下方之程式碼，計算出：

1. 習題一：在一步態週期過程中的 Rg2p、Vg2p、Pcoord_local、Rg2t、Vg2t、Tcoord_local、
Rg2s、Vg2s、Scoord_local、Rg2f、Vg2f、Fcoord_local；
2. 習題二：在一步態週期過程中，右側的膝關節中心與踝關節中心相對於右腳大腿的局部座標系統
的相對運動，並分別命名為 knee_local 和 ankle_local。
'''
Pcoord, Tcoord, Scoord, Fcoord = prepare(vicon_data,start_frame,end_frame)
rRg2p, rVg2p, rPcoord_local = CoordPelvis(Pcoord, 'r')
rRg2t, rVg2t, rTcoord_local = CoordThigh(Tcoord, 'r')
rRg2s, rVg2s, rScoord_local = CoordShank(Scoord, 'r')
rRg2f, rVg2f, rFcoord_local = CoordFoot(Fcoord, 'r')
lRg2p, lVg2p, lPcoord_local = CoordPelvis(Pcoord, 'l')
lRg2t, lVg2t, lTcoord_local = CoordThigh(Tcoord, 'l')
lRg2s, lVg2s, lScoord_local = CoordShank(Scoord, 'l')
lRg2f, lVg2f, lFcoord_local = CoordFoot(Fcoord, 'l')


P_global = {'RLFC':Tcoord["RLFC"],'RMFC':Tcoord["RMFC"],'RlMA':Scoord["RLMA"],'RMMA':Scoord["RMMA"]}
P_local = CoordG2L(rRg2t, rVg2t, P_global)
P_global1 = CoordL2G(rRg2t, rVg2t, P_local)

rknee_local = P_local["RKJC"]
rankle_local = P_local["RAJC"]




'''
下方程式碼將把上方計算之變數存為 .npy 檔，以利批改。
'''


ex1_answer = {'rRg2p':rRg2p, 'rVg2p':rVg2p, 'rPcoord_local':rPcoord_local, 'rRg2t':rRg2t, 'rVg2t':rVg2t, 'rTcoord_local':rTcoord_local, 'rRg2s':rRg2s, 'rVg2s':rVg2s, 'rScoord_local':rScoord_local, 'rRg2f':rRg2f, 'rVg2f':rVg2f, 'rFcoord_local':rFcoord_local, 'lRg2p':lRg2p, 'lVg2p':lVg2p, 'lPcoord_local':lPcoord_local, 'lRg2t':lRg2t, 'lVg2t':lVg2t, 'lTcoord_local':lTcoord_local, 'lRg2s':lRg2s, 'lVg2s':lVg2s, 'lScoord_local':lScoord_local, 'lRg2f':lRg2f, 'lVg2f':lVg2f, 'lFcoord_local':lFcoord_local}
ex2_answer = {'rknee_local':rknee_local, 'rankle_local':rankle_local}
np.save('ex1_answer.npy', ex1_answer)
np.save('ex2_answer.npy', ex2_answer)