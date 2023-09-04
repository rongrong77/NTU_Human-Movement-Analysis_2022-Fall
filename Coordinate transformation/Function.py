import numpy as np
import c3d


def readc3d(file_path):
    # readc3d: 用來讀取c3d檔（輸出為list的版本）
    reader = c3d.Reader(open(file_path, 'rb'))  # load c3d file
    markers = {}
    analog = {}
    marker_name = list(map(lambda name: name.rstrip(), np.ndarray.tolist(reader.point_labels)))
    analog_name = list(map(lambda name: name.rstrip(), np.ndarray.tolist(reader.analog_labels)))
    analog_rate = reader.analog_rate
    frame_rate = reader.point_rate
    for i, j, k in reader.read_frames():  # i=frame, j=markers, k=analog
        if i == 1:
            for n in range(len(marker_name)):
                markers[marker_name[n]] = [j[n].tolist()]
            for n in range(len(analog_name)):
                analog[analog_name[n]] = k[n].tolist()
        else:
            for n in range(len(marker_name)):
                markers[marker_name[n]].append(j[n].tolist())
            for n in range(len(analog_name)):
                analog[analog_name[n]] = analog[analog_name[n]] + k[n].tolist()
    output = markers.copy()
    output.update(analog)
    output_name = marker_name + analog_name
    return output, output_name, frame_rate, analog_rate


def gait_cycle(BW, vicon_data, frame_rate, analog_rate):
    # gait_cycle: 計算步態週期的起始幀和結束幀（三塊力板版本）
    event = []
    for frame in range(len(vicon_data['Fz1'])):
        if -vicon_data['Fz1'][frame] > BW*9.81*0.05:
            event.append(frame)
            break
    for frame in range(len(vicon_data['Fz3'])):
        if -vicon_data['Fz3'][frame] > BW*9.81*0.05:
            event.append(frame)
            break
    return round(min(event)*frame_rate/analog_rate), round(max(event)*frame_rate/analog_rate)


'''
從這裡開始為各位作業一需撰寫的函式
包含習題一的函式 CoordPelvis、CoordThigh、CoordShank、CoordFoot
以及習題二的函式 CoordG2L、CoordL2G
請將下列六個函式的 statement 改寫為作業要求的形式
'''


def CoordPelvis(Pcoord, side):
    # CoordPelvis: 計算骨盆相對global座標系統之旋轉矩陣與位置向量
    Rg2p = 0
    Vg2p = 0
    Pcoord_local = 0
    return Rg2p, Vg2p, Pcoord_local


def CoordThigh(Tcoord, side):
    # CoordThigh: 計算大腿相對global座標系統之旋轉矩陣與位置向量
    Rg2t = 0
    Vg2t = 0
    Tcoord_local = 0
    return Rg2t, Vg2t, Tcoord_local


def CoordShank(Scoord, side):
    # CoordShank: 計算小腿相對global座標系統之旋轉矩陣與位置向量
    Rg2s = 0
    Vg2s = 0
    Scoord_local = 0
    return Rg2s, Vg2s, Scoord_local


def CoordFoot(Fcoord, side):
    # CoordFoot: 計算足部相對global座標系統之旋轉矩陣與位置向量
    Rg2f = 0
    Vg2f = 0
    Fcoord_local = 0
    return Rg2f, Vg2f, Fcoord_local


def CoordG2L(Rg2l, Vg2l, P_global):
    # CoordG2L: 將global座標值轉換至local座標系統對應之座標值
    P_local = 0
    return P_local


def CoordL2G(Rg2l, Vg2l, P_local):
    # CoordL2G: 將local座標值轉換至global座標系統對應之座標值
    P_global = 0
    return P_global