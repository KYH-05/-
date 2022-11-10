#화학세특
#H2
#옥텟의 확장은 고려 안함
#간단한 구조만
#png도 다양하게
#png에 뭐가 뭔지 나타내기
#-------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
#-------------------------------------------------------------------------------
ca=0#central atom(최종결정 중심원자)
cai=[]#central atom index(원자가전자가 4에 가까운 모든 원자의 위치)
mc=[]#molecular components(입력받는 분자구성 요소)
venlist=[]#원자가전자 비교 리스트
elist=[]#전기음성도 비교 리스트
ven={"H":1,"Li":1,"Be":2,"B":3,"C":4,"N":5,"O":6,"F":7,"Na":1,"Mg":2,"Al":3,"Si":4,"P":5,"S":6,"Cl":7,"K":1,"Ca":2,"Ga":3,"Ge":4,"As":5,"Se":6,"Br":7,"Rb":1,"Sr":2,"In":3,"Sn":4,"Sb":5,"Te":6,"I":7,}#valence electron number(원자가전자)
un={"H":1,"Li":1,"Be":2,"B":3,"C":4,"N":3,"O":2,"F":1,"Na":1,"Mg":2,"Al":3,"Si":4,"P":3,"S":2,"Cl":1,"K":1,"Ca":2,"Ga":3,"Ge":4,"As":3,"Se":2,"Br":1,"Rb":1,"Sr":2,"In":3,"Sn":4,"Sb":3,"Te":2,"I":1,}#unpaired electron(홀전자)
e={"Li":0.98,"Be":1.57,"B":2.04,"C":2.55,"N":3.04,"O":3.44,"F":3.98,"Na":0.93,"Mg":1.31,"Al":1.61,"Si":1.9,"P":2.19,"S":2.58,"Cl":3.16,"K":0.82,"Ca":1,"Ga":1.81,"Ge":2.01,"As":2.18,"Se":2.55,"Br":2.96,"Rb":0.82,"Sr":0.95,"In":1.78,"Sn":1.96,"Sb":2.05,"Te":2.1,"I":2.66,}#electronegavity(전기음성도)
#-------------------------------------------------------------------------------
print("분자식의 구성 요소를 하나씩 작성하고 완료되면 X를 쳐주세요.")
print("-"*30)
print("예) H₂O \n->H \n->H \n->O \n->X")
print("-"*30)
while True:
  mc.append(input("->"))
  if "X" in mc:
    break    
mc.remove("X")
#-------------------------------------------------------------------------------
#예외
if len(mc)==2 and mc[0]=="H" and mc[1]=="H":
  print("이원자분자")
  plt.figure(1)
  plt.text(0.435,0.42,"H",fontsize = 60,color='blue')#중심원자
  plt.text(0.16,0.42,"H",fontsize = 60,color='blue')#나머지원자 
  plt.text(0.3,0.45,"―",fontsize = 45,color='blue')
  plt.axis('off') 
  plt.xlim([0, 1])     
  plt.ylim([0, 1])
  plt.figure(2)
  image = img.imread('./이원자분자.png') 
  plt.imshow(image) 
  plt.axis('off') 
  plt.show()  
  input()
#-------------------------------------------------------------------------------
#중심원자 찾기
print(mc)
venlist=[]
elist=[]
for i in range(0,len(mc)):
  if mc[i]=="H":
    venlist.append(100)
  else:
    venlist.append(abs(ven[mc[i]]-4))
cai= [i for i, value in enumerate(venlist) if value == min(venlist)] 
if len(cai)==1:
  ca=mc[cai[0]]
else:
  for i in range(0,len(cai)):
    elist.append(e[mc[cai[i]]])
  ca=mc[cai[elist.index(min(elist))]]
print("중심원자",ca)
#-------------------------------------------------------------------------------
#1.원자가 전자수 합
value=0
for i in range(0,len(mc)):
  value=value+ven[mc[i]]
#-------------------------------------------------------------------------------
#2.기본 골격
mc.remove(ca)
a_sie,b_sie,c_sie,d_sie=0,0,0,0#shared inner electronics
sie=[a_sie,b_sie,c_sie,d_sie]
for i in range(0,len(mc)):
  sie[i]=2
  value=value-2
#-------------------------------------------------------------------------------
#3.바같원자 옥텟규칙
a_nsoe,b_nsoe,c_nsoe,d_nsoe=0,0,0,0#non shared outside electronics
nsoe=[a_nsoe,b_nsoe,c_nsoe,d_nsoe]
for i in range(0,len(mc)):
  if mc[i]=="H" or mc[i]=="Li" or mc[i]=="Na" or mc[i]=="K" or mc[i]=="Rb":
    nsoe[i]=0
  else:
    nsoe[i]=6
    value=value-6
#-------------------------------------------------------------------------------
#4.중심원자 나머지 채우기
nsie=value#non shard inner electronics
#-------------------------------------------------------------------------------
#5.중심원자 옥텟규칙-Be,B가 옥텟 예외상황도 고려
ie=nsie+sie[0]+sie[1]+sie[2]+sie[3]#inner electronics
if ven[ca]==2:#2족
  if ie==4:
    print("종료")
elif ven[ca]==3:#13족
  if ie==6:
    print("종료")
else:#정상
  if ie==8:#단일결합
    print("종료")
  else:#다중결합
    ien=8-ie#inner electronics need
    egi=[]#electronics giver index
    for i in range(0,len(nsoe)):
      if nsoe[i]>0:
        egi.append(i)
    if len(egi)==1:#줄수있는 원자가 1개
      nsoe[egi[0]]=nsoe[egi[0]]-ien
      sie[egi[0]]=sie[egi[0]]+ien
      print("종료")
    else:#줄수 있는 원자가 2개 이상-아마 중심원소가 원자가전자4인 경우만
      #nsoe[0]=nsoe[0]-2
      #sie[0]=sie[0]+2
      #nsoe[1]=nsoe[1]-2
      #sie[1]=sie[1]+2
      #print("종료")
      fc=[[],[],[]]#Formal charge형식전하 리스트
      fc_c=[]#형식전하 비교 리스트
      if len(mc)==2:
        #31
        sie[0]=sie[0]+4
        nsoe[0]=nsoe[0]-4
        fc[0].append(ven[mc[0]]-nsoe[0]-(1/2)*sie[0])
        fc[0].append(ven[ca]-nsie-(1/2)*(sie[0]+sie[1]+sie[2]+sie[3]))
        fc[0].append(ven[mc[1]]-nsoe[1]-(1/2)*sie[1])
        sie[0]=sie[0]-4
        nsoe[0]=nsoe[0]+4
        #13
        sie[1]=sie[1]+4
        nsoe[1]=nsoe[1]-4
        fc[1].append(ven[mc[0]]-nsoe[0]-(1/2)*sie[0])
        fc[1].append(ven[ca]-nsie-(1/2)*(sie[0]+sie[1]+sie[2]+sie[3]))
        fc[1].append(ven[mc[1]]-nsoe[1]-(1/2)*sie[1])
        sie[1]=sie[1]-4
        nsoe[1]=nsoe[1]+4  
        #22
        nsoe[0]=nsoe[0]-2
        sie[0]=sie[0]+2
        nsoe[1]=nsoe[1]-2
        sie[1]=sie[1]+2          
        fc[2].append(ven[mc[0]]-nsoe[0]-(1/2)*sie[0])
        fc[2].append(ven[ca]-nsie-(1/2)*(sie[0]+sie[1]+sie[2]+sie[3]))
        fc[2].append(ven[mc[1]]-nsoe[1]-(1/2)*sie[1])
        nsoe[0]=nsoe[0]+2
        sie[0]=sie[0]-2
        nsoe[1]=nsoe[1]+2
        sie[1]=sie[1]-2    
        #계산
        for i in range(0,3):
          fc_c.append((fc[i][0])**2+(fc[i][1])**2+(fc[i][2])**2)
        fc_c_r=fc_c.copy()
        fc_c.sort()
        if fc_c[0]< fc_c[1]:
          a=fc_c_r.index(fc_c[0])
          if a==0: #31
            sie[0]=sie[0]+4
            nsoe[0]=nsoe[0]-4
            print("종료")
          elif a==1:#13
            sie[1]=sie[1]+4
            nsoe[1]=nsoe[1]-4
            print("종료")
          elif a==2:#22
            nsoe[0]=nsoe[0]-2
            sie[0]=sie[0]+2
            nsoe[1]=nsoe[1]-2
            sie[1]=sie[1]+2    
            print("종료")
        elif fc_c[0]==fc_c[1]:#전기음성도 계산 필요
          if e[mc[0]]>e[mc[1]]:#전기음성도가 큰 쪽이 단일결합
            nsoe[1]=nsoe[1]-4
            sie[1]=sie[1]+4
            print("종료")
          if e[mc[0]]<e[mc[1]]:
            nsoe[0]=nsoe[0]-4
            sie[0]=sie[0]+4     
            print("종료")
#-------------------------------------------------------------------------------          
#루이스 전자점식
#원자----------------- 
plt.figure(1)
real_lenmc=len(mc)
for i in range(0,4-len(mc)):
  mc.append("")
plt.text(0.435,0.42,ca,fontsize = 60,color='blue')#중심원자
plt.text(0.16,0.42,mc[0],fontsize = 60,color='blue')#나머지원자
plt.text(0.71,0.42,mc[1],fontsize = 60,color='blue')
plt.text(0.435,0.760,mc[2],fontsize = 60,color='blue')
plt.text(0.435,0.085,mc[3],fontsize = 60,color='blue')
#점----------------- 
if real_lenmc>=1:
  if mc[0]!="H":
    if nsoe[0]==2:
      plt.text(0.1,0.5,"•",fontsize = 30,color='blue')
      plt.text(0.1,0.42,"•",fontsize = 30,color='blue')      
    elif nsoe[0]==4:
      plt.text(0.1,0.5,"•",fontsize = 30,color='blue')
      plt.text(0.1,0.42,"•",fontsize = 30,color='blue')  
      plt.text(0.165,0.59,"•",fontsize = 30,color='blue')
      plt.text(0.23,0.59,"•",fontsize = 30,color='blue')      
    elif nsoe[0]==6:
      plt.text(0.165,0.59,"•",fontsize = 30,color='blue')
      plt.text(0.23,0.59,"•",fontsize = 30,color='blue')
      plt.text(0.1,0.5,"•",fontsize = 30,color='blue')
      plt.text(0.1,0.42,"•",fontsize = 30,color='blue')
      plt.text(0.165,0.33,"•",fontsize = 30,color='blue')
      plt.text(0.23,0.33,"•",fontsize = 30,color='blue')
if real_lenmc>=2:
  if mc[0]!="H":
    if nsoe[1]==2:
      plt.text(0.83,0.5,"•",fontsize = 30,color='blue')
      plt.text(0.83,0.42,"•",fontsize = 30,color='blue')      
    elif nsoe[1]==4:
      plt.text(0.83,0.5,"•",fontsize = 30,color='blue')
      plt.text(0.83,0.42,"•",fontsize = 30,color='blue')  
      plt.text(0.71,0.59,"•",fontsize = 30,color='blue')
      plt.text(0.782,0.59,"•",fontsize = 30,color='blue')      
    elif nsoe[1]==6:
      plt.text(0.71,0.59,"•",fontsize = 30,color='blue')
      plt.text(0.782,0.59,"•",fontsize = 30,color='blue')
      plt.text(0.83,0.5,"•",fontsize = 30,color='blue')
      plt.text(0.83,0.42,"•",fontsize = 30,color='blue')
      plt.text(0.71,0.33,"•",fontsize = 30,color='blue')
      plt.text(0.782,0.33,"•",fontsize = 30,color='blue')
if real_lenmc>=3:
  if mc[0]!="H":
    if nsoe[2]==2:
      plt.text(0.435,0.93,"•",fontsize = 30,color='blue')
      plt.text(0.507,0.93,"•",fontsize = 30,color='blue')      
    elif nsoe[2]==4:
      plt.text(0.435,0.93,"•",fontsize = 30,color='blue')
      plt.text(0.507,0.93,"•",fontsize = 30,color='blue')
      plt.text(0.56,0.85,"•",fontsize = 30,color='blue')
      plt.text(0.56,0.77,"•",fontsize = 30,color='blue')      
    elif nsoe[2]==6:
      plt.text(0.435,0.93,"•",fontsize = 30,color='blue')
      plt.text(0.507,0.93,"•",fontsize = 30,color='blue')
      plt.text(0.38,0.85,"•",fontsize = 30,color='blue')
      plt.text(0.38,0.77,"•",fontsize = 30,color='blue')
      plt.text(0.56,0.85,"•",fontsize = 30,color='blue')
      plt.text(0.56,0.77,"•",fontsize = 30,color='blue')
if real_lenmc>=4:
  if mc[0]!="H":
    if nsoe[3]==2:
      plt.text(0.435,0.005,"•",fontsize = 30,color='blue')
      plt.text(0.507,0.005,"•",fontsize = 30,color='blue')        
    elif nsoe[3]==4:
      plt.text(0.435,0.005,"•",fontsize = 30,color='blue')
      plt.text(0.507,0.005,"•",fontsize = 30,color='blue')    
      plt.text(0.38,0.17,"•",fontsize = 30,color='blue')
      plt.text(0.38,0.085,"•",fontsize = 30,color='blue')      
    elif nsoe[3]==6:
      plt.text(0.435,0.005,"•",fontsize = 30,color='blue')
      plt.text(0.507,0.005,"•",fontsize = 30,color='blue')  
      plt.text(0.38,0.17,"•",fontsize = 30,color='blue')
      plt.text(0.38,0.085,"•",fontsize = 30,color='blue')
      plt.text(0.56,0.17,"•",fontsize = 30,color='blue')
      plt.text(0.56,0.085,"•",fontsize = 30,color='blue')
#선----------------- 
if real_lenmc>=1:
  if sie[0]==2:
    plt.text(0.3,0.45,"―",fontsize = 45,color='blue')
  elif sie[0]==4:
    plt.text(0.3,0.41,"―",fontsize = 45,color='blue')
    plt.text(0.3,0.49,"―",fontsize = 45,color='blue')
  elif sie[0]==6:
    plt.text(0.3,0.49,"―",fontsize = 45,color='blue')
    plt.text(0.3,0.45,"―",fontsize = 45,color='blue')
    plt.text(0.3,0.41,"―",fontsize = 45,color='blue')
if real_lenmc>=2:
  if sie[1]==2:
    plt.text(0.575,0.45,"―",fontsize = 45,color='blue')
  elif sie[1]==4:
    plt.text(0.575,0.41,"―",fontsize = 45,color='blue')
    plt.text(0.575,0.49,"―",fontsize = 45,color='blue')
  elif sie[1]==6:
    plt.text(0.575,0.49,"―",fontsize = 45,color='blue')
    plt.text(0.575,0.45,"―",fontsize = 45,color='blue')
    plt.text(0.575,0.41,"―",fontsize = 45,color='blue')
if real_lenmc>=3:
  if sie[2]==2:
    plt.text(0.48,0.63,"|",fontsize = 42,color='blue')
  elif sie[2]==4:
    plt.text(0.44,0.63,"|",fontsize = 42,color='blue')
    plt.text(0.52,0.63,"|",fontsize = 42,color='blue')    
  elif sie[2]==6:
    plt.text(0.44,0.63,"|",fontsize = 42,color='blue')
    plt.text(0.48,0.63,"|",fontsize = 42,color='blue')
    plt.text(0.52,0.63,"|",fontsize = 42,color='blue')
if real_lenmc>=4:
  if sie[3]==2:
    plt.text(0.48,0.29,"|",fontsize = 42,color='blue')
  elif sie[3]==4:
    plt.text(0.44,0.29,"|",fontsize = 42,color='blue')
    plt.text(0.52,0.29,"|",fontsize = 42,color='blue')    
  elif sie[3]==6:  
    plt.text(0.44,0.29,"|",fontsize = 42,color='blue')
    plt.text(0.48,0.29,"|",fontsize = 42,color='blue')
    plt.text(0.52,0.29,"|",fontsize = 42,color='blue')
#내부점----------------- 
if nsie==2:
  plt.text(0.435,0.33,"•",fontsize = 30,color='blue')
  plt.text(0.507,0.33,"•",fontsize = 30,color='blue')
elif nsie==4:
  plt.text(0.435,0.33,"•",fontsize = 30,color='blue')
  plt.text(0.507,0.33,"•",fontsize = 30,color='blue') 
  plt.text(0.435,0.59,"•",fontsize = 30,color='blue')
  plt.text(0.507,0.59,"•",fontsize = 30,color='blue')  
elif nsie==6:
  plt.text(0.435,0.33,"•",fontsize = 30,color='blue')
  plt.text(0.507,0.33,"•",fontsize = 30,color='blue') 
  plt.text(0.435,0.59,"•",fontsize = 30,color='blue')###################################### HCl오류 여기에 있을듯
  plt.text(0.507,0.59,"•",fontsize = 30,color='blue') 
  plt.text(0.58,0.53,"•",fontsize = 30,color='blue')
  plt.text(0.58,0.4,"•",fontsize = 30,color='blue')      
#기타----------------- 
plt.axis('off') 
plt.xlim([0, 1])     
plt.ylim([0, 1])
#-------------------------------------------------------------------------------
#형태 출력,입체 모델 구현
plt.figure(2)
if real_lenmc==1:
  image = img.imread('./이원자분자.png') 
  print("이원자분자")
elif real_lenmc==2:
  if nsie==0:
    image = img.imread('./선형.png') 
    print("직선형")
  else:
    image = img.imread('./굽은형.png') 
    print("굽은형")
elif real_lenmc==3:  
  if nsie==0:
    image = img.imread('./삼각형.png') 
    print("평면 삼각형")
  else:
    image = img.imread('./삼각뿔.png') 
    print("삼각뿔")
else:
  image = img.imread('./사면체.png') 
  print("사면체")
plt.imshow(image) 
plt.axis('off') 
plt.show()

