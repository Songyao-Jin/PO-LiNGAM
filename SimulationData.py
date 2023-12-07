import numpy as np
import pandas as pd

def ToBij():
     ten = np.random.randint(0,2)
     s  =np.random.random()
     while abs(s) <0.5 and ten ==0:
        s  =np.random.random()
     result = ten+s
     if np.random.randint(0,10)>5:
        result = -1*result
     return round(result,3)


def SelectPdf(Num,k=3):
    noise = pow(np.random.exponential(scale=1, size=Num),2)
    return noise


def GaussianPdf(Num):
    noise = np.random.normal(0, 1, Num)
    return noise


def Toa():
    return 0.2




def CaseS1(Num=5000):
    x1=SelectPdf(Num)
    x2=SelectPdf(Num)+ToBij()*x1
    x3=SelectPdf(Num)+ToBij()*x1+ToBij()*x2
    x4=SelectPdf(Num)+ToBij()*x1+ToBij()*x3
    x5=SelectPdf(Num)+ToBij()*x1+ToBij()*x2
    x6=SelectPdf(Num)+ToBij()*x2+ToBij()*x3
    
    
    data = pd.DataFrame(np.array([x1,x2,x3,x4,x5,x6]).T,columns=['x1','x2','x3','x4','x5','x6'])

    data = (data-data.mean())/data.std()
    #data = data-data.mean()
    return data  




def CaseS2(Num=5000):
    x1=SelectPdf(Num)
    L1=SelectPdf(Num)+ToBij()*x1
    L2=SelectPdf(Num)+ToBij()*x1+ToBij()*L1
    x2=SelectPdf(Num)+ToBij()*L1
    x3=SelectPdf(Num)+ToBij()*L1
    x4=SelectPdf(Num)+ToBij()*L1
    
    x5=SelectPdf(Num)+ToBij()*L2
    x6=SelectPdf(Num)+ToBij()*L2      
    x7=SelectPdf(Num)+ToBij()*L2        
    x8=SelectPdf(Num)+ToBij()*L2+ToBij()*x7
    
    data = pd.DataFrame(np.array([x1,x2,x3,x4,x5,x6,x7,x8]).T,columns=['x1','x2','x3','x4','x5','x6','x7','x8'])

    data = (data-data.mean())/data.std()
    #data = data-data.mean()
    return data  



def CaseS3(Num=5000):
    L1=SelectPdf(Num)
    L2=SelectPdf(Num)+ToBij()*L1
    L3=SelectPdf(Num)+ToBij()*L1
    L4=SelectPdf(Num)+ToBij()*L1
    
    x1=SelectPdf(Num)+ToBij()*L2
    x2=SelectPdf(Num)+ToBij()*L2
    
    x3=SelectPdf(Num)+ToBij()*L2+ToBij()*L3  
    x4=SelectPdf(Num)+ToBij()*L3 
    x5=SelectPdf(Num)+ToBij()*L3 
    x6=SelectPdf(Num)+ToBij()*L4       
    x7=SelectPdf(Num)+ToBij()*L4
    x8=SelectPdf(Num)+ToBij()*L4       
    
    data = pd.DataFrame(np.array([x1,x2,x3,x4,x5,x6,x7,x8]).T,columns=['x1','x2','x3','x4','x5','x6','x7','x8'])

    data = (data-data.mean())/data.std()
    #data = data-data.mean()
    return data  




def main():
    data = CaseS1()
    print(data['x1'])

if __name__ == '__main__':
    main()
