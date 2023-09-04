import numpy as np
import c3d
import math as m
import sympy as sp


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
    
    Vg2p=Pcoord[0]
    Pcoord_local=0
    Rg2p=np.empty([169,3,3])
    
    if side =='r':
    
        for i in range(Pcoord.shape[1]):
            zp=(Pcoord[0][i]-Pcoord[1][i])/np.linalg.norm(Pcoord[0][i]-Pcoord[1][i]) 
            yp=np.cross(zp,(Pcoord[0][i]-Pcoord[2][i]))/np.linalg.norm(np.cross(zp,(Pcoord[0][i]-Pcoord[2][i])))
            xp= np.cross(yp,zp)
            
    
            Rg2p[i,:,0]=xp
            Rg2p[i,:,1]=yp
            Rg2p[i,:,2]=zp
        
            Pcoord_local =np.empty([3,169,3])
        
            for j in range(Pcoord.shape[1]):
                Pcoord_local[0,j,:] = np.dot(np.transpose(Rg2p[j,:,:]),Pcoord[0,j,:]-Vg2p[j,:])
                Pcoord_local[1,j,:] = np.dot(np.transpose(Rg2p[j,:,:]),Pcoord[1,j,:]-Vg2p[j,:])
                Pcoord_local[2,j,:] = np.dot(np.transpose(Rg2p[j,:,:]),Pcoord[2,j,:]-Vg2p[j,:])

    elif side=='l':
        for i in range(Pcoord.shape[1]):
            zp=-(Pcoord[0][i]-Pcoord[1][i])/np.linalg.norm(Pcoord[0][i]-Pcoord[1][i]) 
            yp=np.cross(zp,(Pcoord[0][i]-Pcoord[2][i]))/np.linalg.norm(np.cross(zp,(Pcoord[0][i]-Pcoord[2][i])))
            xp= np.cross(yp,zp)
          
            Rg2p[i,:,0]=xp
            Rg2p[i,:,1]=yp
            Rg2p[i,:,2]=zp 
         
            Pcoord_local =np.empty([3,169,3])
           
            for j in range(Pcoord.shape[1]):
                Pcoord_local[0,j,:] = np.dot(np.transpose(Rg2p[j,:,:]),Pcoord[0,j,:]-Vg2p[j,:])
                Pcoord_local[1,j,:] = np.dot(np.transpose(Rg2p[j,:,:]),Pcoord[1,j,:]-Vg2p[j,:])
                Pcoord_local[2,j,:] = np.dot(np.transpose(Rg2p[j,:,:]),Pcoord[2,j,:]-Vg2p[j,:])
                
    return Rg2p, Vg2p, Pcoord_local



def CoordThigh(Tcoord, side):
    # CoordThigh: 計算大腿相對global座標系統之旋轉矩陣與位置向量
    
    Vg2t=Tcoord[0]
    Rg2t =np.empty([169,3,3])
    
    if side =='r':
    
        for i in range(Tcoord.shape[1]):
            zrt=(Tcoord[1][i]-Tcoord[2][i])/np.linalg.norm(Tcoord[1][i]-Tcoord[2][i])
            xrt=np.cross((Tcoord[0][i]-Tcoord[1][i]),zrt)/np.linalg.norm(np.cross((Tcoord[0][i]-Tcoord[1][i]),zrt))
            yrt= np.cross(zrt,xrt)
    
            Rg2t[i,:,0]=xrt
            Rg2t[i,:,1]=yrt
            Rg2t[i,:,2]=zrt
            
            Tcoord_local =np.empty([3,169,3])
            
            for j in range(169):
                Tcoord_local[0,j,:] = np.dot(np.transpose(Rg2t[j,:,:]),Tcoord[0,j,:]-Vg2t[j,:])
                Tcoord_local[1,j,:] = np.dot(np.transpose(Rg2t[j,:,:]),Tcoord[1,j,:]-Vg2t[j,:])
                Tcoord_local[2,j,:] = np.dot(np.transpose(Rg2t[j,:,:]),Tcoord[2,j,:]-Vg2t[j,:])
               

                
    elif side=='l':
        
        for i in range(Tcoord.shape[1]):
            zrt=(Tcoord[1][i]-Tcoord[2][i])/np.linalg.norm(Tcoord[1][i]-Tcoord[2][i])
            xrt=np.cross((Tcoord[0][i]-Tcoord[1][i]),zrt)/np.linalg.norm(np.cross((Tcoord[0][i]-Tcoord[1][i]),zrt))
            yrt= np.cross(zrt,xrt)
    
          
            Rg2t[i,:,0]=xrt
            Rg2t[i,:,1]=yrt
            Rg2t[i,:,2]=zrt
         
            Tcoord_local =np.empty([3,169,3])
           
            for j in range(Tcoord.shape[1]):
                Tcoord_local[0,j,:] = np.dot(np.transpose(Rg2t[j,:,:]),Tcoord[0,j,:]-Vg2t[j,:])
                Tcoord_local[1,j,:] = np.dot(np.transpose(Rg2t[j,:,:]),Tcoord[1,j,:]-Vg2t[j,:])
                Tcoord_local[2,j,:] = np.dot(np.transpose(Rg2t[j,:,:]),Tcoord[2,j,:]-Vg2t[j,:])
    
    return Rg2t, Vg2t, Tcoord_local


def CoordShank(Scoord, side):
    # CoordShank: 計算小腿相對global座標系統之旋轉矩陣與位置向量

    Vg2s=Scoord[0]
    Scoord_local=0
    Rg2s =np.empty([169,3,3])
    
    if side =='r':
    
        for i in range(Scoord.shape[1]):
            xrs=np.divide(np.cross(Scoord[1][i]-Scoord[3][i],Scoord[2][i]-Scoord[3][i]),np.linalg.norm(np.cross(Scoord[1][i]-Scoord[3][i],Scoord[2][i]-Scoord[3][i])))
            zrs=np.divide(np.cross(xrs,(Scoord[0][i]-(Scoord[3][i]+Scoord[2][i])/2)),np.linalg.norm(np.cross(xrs,Scoord[0][i]-(Scoord[3][i]+Scoord[2][i])/2)))
            yrs= np.cross(zrs,xrs)
            
            Rg2s[i,:,0]=xrs
            Rg2s[i,:,1]=yrs
            Rg2s[i,:,2]=zrs
            
            Scoord_local =np.empty([4,169,3])

            for j in range(Scoord.shape[1]):
                Scoord_local[0,j,:] = np.dot(np.transpose(Rg2s[j,:,:]),Scoord[0,j,:]-Vg2s[j,:])
                Scoord_local[1,j,:] = np.dot(np.transpose(Rg2s[j,:,:]),Scoord[1,j,:]-Vg2s[j,:])
                Scoord_local[2,j,:] = np.dot(np.transpose(Rg2s[j,:,:]),Scoord[2,j,:]-Vg2s[j,:])
                Scoord_local[3,j,:] = np.dot(np.transpose(Rg2s[j,:,:]),Scoord[3,j,:]-Vg2s[j,:])
                 
    elif side=='l':
        
        for i in range(Scoord.shape[1]):
            xrs=np.divide(np.cross(Scoord[1][i]-Scoord[3][i],Scoord[2][i]-Scoord[3][i]),np.linalg.norm(np.cross(Scoord[1][i]-Scoord[3][i],Scoord[2][i]-Scoord[3][i])))
            zrs=np.divide(np.cross(xrs,(Scoord[0][i]-0.5*(Scoord[3][i]+Scoord[2][i]))),np.linalg.norm(np.cross(xrs,Scoord[0][i]-0.5*(Scoord[3][i]+Scoord[2][i]))))
            yrs= np.cross(zrs,xrs)
            
            Rg2s[i,:,0]=xrs
            Rg2s[i,:,1]=yrs
            Rg2s[i,:,2]=zrs
            
            Scoord_local =np.empty([4,169,3])

            for j in range(Scoord.shape[1]):
                Scoord_local[0,j,:] = np.dot(np.transpose(Rg2s[j,:,:]),Scoord[0,j,:]-Vg2s[j,:])
                Scoord_local[1,j,:] = np.dot(np.transpose(Rg2s[j,:,:]),Scoord[1,j,:]-Vg2s[j,:])
                Scoord_local[2,j,:] = np.dot(np.transpose(Rg2s[j,:,:]),Scoord[2,j,:]-Vg2s[j,:])
                Scoord_local[3,j,:] = np.dot(np.transpose(Rg2s[j,:,:]),Scoord[3,j,:]-Vg2s[j,:])
            

    return Rg2s, Vg2s, Scoord_local

def CoordFoot(Fcoord, side):
    # CoordFoot: 計算足部相對global座標系統之旋轉矩陣與位置向量
    
    Vg2f=Fcoord[0]
    Fcoord_local=0
    Rg2f =np.empty([169,3,3])
    
    if side =='r':
    
        for i in range(Fcoord.shape[1]):
            xrf =((Fcoord[1][i]+Fcoord[2][i])/2-Fcoord[0][i])/np.linalg.norm((Fcoord[1][i]+Fcoord[2][i])/2-Fcoord[0][i])
            yrf =np.cross(xrf,(Fcoord[1][i]-Fcoord[2][i]))/np.linalg.norm(np.cross(xrf,(Fcoord[1][i]-Fcoord[2][i])))
            zrf =np.cross(xrf,yrf)
    
            Rg2f[:,:,0]=xrf
            Rg2f[:,:,1]=yrf
            Rg2f[:,:,2]=zrf
            
            Fcoord_local =np.empty([3,169,3])

            for j in range(Fcoord.shape[1]):
                Fcoord_local[0,j,:] = np.dot(np.transpose(Rg2f[j,:,:]),Fcoord[0,j,:]-Vg2f[j,:])
                Fcoord_local[1,j,:] = np.dot(np.transpose(Rg2f[j,:,:]),Fcoord[1,j,:]-Vg2f[j,:])

    elif side =='l':
    
        for i in range(Fcoord.shape[1]):
            xrf =((Fcoord[1][i]+Fcoord[2][i])/2-Fcoord[0][i])/np.linalg.norm((Fcoord[1][i]+Fcoord[2][i])/2-Fcoord[0][i])
            yrf =np.cross(xrf,(Fcoord[1][i]-Fcoord[2][i]))/np.linalg.norm(np.cross(xrf,(Fcoord[1][i]-Fcoord[2][i])))
            zrf =np.cross(xrf,yrf)
    
            Rg2f[:,:,0]=xrf
            Rg2f[:,:,1]=yrf
            Rg2f[:,:,2]=zrf
            
            Fcoord_local =np.empty([3,169,3])

            for j in range(Fcoord.shape[1]):
                 Fcoord_local[0,j,:] = np.dot(np.transpose(Rg2f[j,:,:]),Fcoord[0,j,:]-Vg2f[j,:])
                 Fcoord_local[1,j,:] = np.dot(np.transpose(Rg2f[j,:,:]),Fcoord[1,j,:]-Vg2f[j,:])
                 Fcoord_local[2,j,:] = np.dot(np.transpose(Rg2f[j,:,:]),Fcoord[2,j,:]-Vg2f[j,:])
     
    return Rg2f, Vg2f, Fcoord_local


def CoordL2G(Rg2t, Vg2t, P_local):
    # CoordL2G: 將local座標值轉換至global座標系統對應之座標值
    
    P_global=np.empty([3,169,3])
    
    for i in range(P_local.shape[1]):
        P_global[0,i,:] = np.dot(Rg2t[i,:,:],P_local[0,i,:])+Vg2t[i,:]
        P_global[1,i,:] = np.dot(Rg2t[i,:,:],P_local[1,i,:])+Vg2t[i,:]
        P_global[2,i,:] = np.dot(Rg2t[i,:,:],P_local[2,i,:])+Vg2t[i,:]

    return P_global

def Rot2Ang(Rot, sequence):
    theta = np.empty([169,3])
    
    for i in range(169):
        if sequence =='zxy':
            theta[i,0]=m.degrees(m.atan(-Rot[i,0,1]/Rot[i,1,1]))
            theta[i,2]=m.degrees(m.atan(-Rot[i,2,0]/Rot[i,2,2]))
            if m.cos(theta[i,1])>0:
                theta[i,1]=m.degrees(m.asin(Rot[i,2,1]))
            else:
                if m.asin(m.sin(theta[i,1]))>0:
                    theta[i,1]=180-m.degrees(m.asin(Rot[i,2,1]))
                else:
                    theta[i,1]=-180-m.degrees(m.asin(Rot[i,2,1]))
        elif sequence =='yxz':
              theta[i,0]=m.degrees(m.atan(Rot[i,0,2]/Rot[i,2,2]))
              theta[i,1]=m.degrees(m.asin(-Rot[i,1,2]))
              theta[i,2]=m.degrees(m.atan(Rot[i,1,0]/Rot[i,1,1]))
            
    return theta 
def cal_deriv(x, y):
    diff_x = []                       # 用来存储x列表中的两数之差
    for i, j in zip(x[:,0,:], x[1::]):  
        diff_x.append(j - i)
 
    diff_y = []                       # 用来存储y列表中的两数之差
    for i, j in zip(y[0::], y[1::]):
        diff_y.append(j - i)  
        
    slopes = []                       # 用来存储斜率
    for i in range(len(diff_y)):
        slopes.append(diff_y[i] / diff_x[i])
        
    deriv = []                        # 用来存储一阶导数
    for i, j in zip(slopes[0::], slopes[1::]):        
        deriv.append((0.5 * (i + j))) # 根据离散点导数的定义，计算并存储结果
    deriv.insert(0, slopes[0])        # (左)端点的导数即为与其最近点的斜率
    deriv.append(slopes[-1])          # (右)端点的导数即为与其最近点的斜率
 
    for i in deriv:                   # 打印结果，方便检查，调用时也可注释掉
        print(i)
        
    return deriv   
def cal_2nd_deriv(x,y):
    return cal_deriv(x,cal_deriv(x,y))
    

def Derivative(x,y,n):
    if n==1:
        dyi=cal_deriv(x,y)
    if n==2:
        dyi=cal_2nd_deriv(x,y)
    return dyi
def Ang2LocalAngular(theta, seq, smprate):
    
    theta = np.empty([169,3])
    AngVel= np.empty([169,3])
    
    if seq =='zxy' :
        for i in range(169):
            R = np.array([[sp.cos(theta[i,2]),0,-sp.cos(theta[i,1])*sp.sin(theta[i,2])],
                          [0, 1,sp.sin(theta[i,1])],
                          [sp.sin(theta[i,2]),0,sp.cos(theta[i,1])*sp.cos(theta[i,2])]])
            t =np.array([[Derivative(theta[i,1], i, 1)],
                        [Derivative(theta[i,2], i, 1)],
                        [Derivative(theta[i,0], i, 1)]])
                         
            AngVel= np.dot(R[i,:,:],t[i,:])
        
  
    return AngVel
'''
def Rot2LocalAngular(Rg2l,smprate):
    
    return AngVel, AngAcc
'''