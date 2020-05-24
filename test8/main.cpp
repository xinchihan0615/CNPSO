
//
//  main.cpp
//  test5
//
//  Created by hanxinchi on 18/5/24.
//  Copyright © 2018年 hanxinchi. All rights reserved.
//

//
//  main.cpp
//  Test(3)
//
//  Created by hanxinchi on 18/5/18.
//  Copyright © 2018年 hanxinchi. All rights reserved.
//

//
//  main.cpp
//  test
//
//  Created by hanxinchi on 18/5/17.
//  Copyright © 2018年 hanxinchi. All rights reserved.
//

//1.chu shi hua partical  2.partical gen xin 3.pbest 4.partical shu xin


#include <time.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include "cec2013.cpp"
#include "cfunction.cpp"

//有关数据集的条件
#define MAX_N 30          //为数组开辟足够的存储空间
int funIndex=9;
int funDim;
int maxFes=400000;
double rhO ;   //数据集要求精度
double maxFitness ;  //数据集提供的最优值
double UB[MAX_N];
double LB[MAX_N];


//有关算法运行的条件
#define initnumOfNiche 10
#define innerRunTime 10
#define popNum 50
#define TimeToSuspendSet 50  //暂停次数阈值
#define accuracy 0.00001  //算法更新过程中的默认精度
int MAX_NICHE_NUM = 10;
double disToSwitch;  //pso与clpso切换的阈值
int totalTime = 0;  //评估次数
int initsizeOfNiche = popNum/initnumOfNiche;

std::default_random_engine e;

using namespace std;
#define PI 3.1415926

//根据测试集要求，设置有关测试的参数
void initVal(int funIndex){
    if((funIndex>=1&&funIndex<=4)||(funIndex>=10&&funIndex<=20))
        rhO = 0.01;
    if(funIndex==5||funIndex==6||funIndex==8)
        rhO = 0.5;
    if(funIndex==7||funIndex==9)
        rhO = 0.2;
    if(funIndex==1){
        for(int i=0;i<funDim;++i){
            LB[i]=0;
            UB[i]=30;
        }
    }
    if(funIndex==2||funIndex==3||funIndex==10){
        for(int i=0;i<funDim;++i){
            LB[i]=0;
            UB[i]=1;
        }
    }
    if(funIndex==4){
        for(int i=0;i<funDim;++i){
            LB[i]=-6;
            UB[i]=6;
        }
    }
    if(funIndex==5){
        LB[0]=-1.9;
        LB[1]=-1.1;
        UB[0]=1.9;
        UB[1]=1.1;
    }
    if(funIndex==6){
        for(int i=0;i<funDim;++i){
            LB[i]=-10;
            UB[i]=10;
        }
    }
    if(funIndex==7){
        for(int i=0;i<funDim;++i){
            LB[i]=0.25;
            UB[i]=10;
        }
    }
    if(funIndex==8){
        for(int i=0;i<funDim;++i){
            LB[i]=-10;
            UB[i]=10;
        }
    }
    if(funIndex==9){
        for(int i=0;i<funDim;++i){
            LB[i]=0.25;
            UB[i]=10;
        }
    }
    if(funIndex>=11){
        for(int i=0;i<funDim;++i){
            LB[i]=-5;
            UB[i]=5;
        }
    }
    if(funIndex==1||funIndex==4)
        maxFitness=200.0;
    if(funIndex==2||funIndex==3||funIndex==7||funIndex==9)
        maxFitness=1.0;
    if(funIndex==5)
        maxFitness=1.03163;
    if(funIndex==6)
        maxFitness=186.7309088310240;
    if(funIndex==8)
        maxFitness=2709.0935;
    if(funIndex==9)
        maxFitness=1.0;
    if(funIndex==10)
        maxFitness=-2;
    if(funIndex==1||funIndex==2||funIndex==3)
        funDim = 1;
    if(funIndex>=4&&funIndex<=7)
        funDim = 2;
    if(funIndex==8||funIndex==9)
        funDim = 3;
    if(funIndex>=10&&funIndex<=13)
        funDim = 2;
    if(funIndex==14||funIndex==15)
        funDim = 3;
    if(funIndex==16||funIndex==17)
        funDim = 5;
    if(funIndex==18||funIndex==19)
        funDim = 10;
    if(funIndex==20)
        funDim = 20;
    
}

void setDisToSwitch()
{
    double disOfDiagonal=0;
    for(int d=0;d<funDim;++d)
    {
        disOfDiagonal+=(UB[d]-LB[d])*(UB[d]-LB[d]);
    }
    disOfDiagonal = sqrt(disOfDiagonal*funDim);
    disToSwitch = disOfDiagonal/50.0;
}

//1. 随机数生成
double r8_uniform_ab(double low, double high) {
    int num = arc4random_uniform(10001);
    return low+num%10000/10000.0*(high - low);
}

//一开始辅粒子随机程度
double rand_back(double i,int d)
{
    double r;
    r=-(UB[d]-LB[d])/15.0+(UB[d]-LB[d])/15.0*2*r8_uniform_ab(0, 1);
    r+=i;
    return r;
}


double rand_back2(double i,int d)
{
    double r;
    r=-(UB[d]-LB[d])/2.0+(UB[d]-LB[d])/2.0*2*r8_uniform_ab(0, 1);
    r+=i;
    r = LB[d]+r8_uniform_ab(0, 1)*(UB[d]-LB[d]);
    //r = r + r8_uniform_ab(0, 1)*(i-r);
    return r;
}




//2. Z值更新
void circleMap(double &z) {
    double b = 11.0;
    if ((z + 0.5 - (b / 2 * PI) * sin(2.0 * PI * z)) - (int)(z + 0.5 - (b / 2 * PI) * sin(2.0 * PI * z))>0) {
        z = (z + 0.5 - (b / 2 * PI) * sin(2.0 * PI * z)) - (int)(z + 0.5 - (b / 2 * PI) * sin(2.0 * PI * z));
    }
    else {
        z = 1+((z + 0.5 - (b / 2 * PI) * sin(2.0 * PI * z)) - (int)(z + 0.5 - (b / 2 * PI) * sin(2.0 * PI * z))) ;
    }
}


//3. 计算函数值
double funTest(int dim, double*x, int EvalCount = -1)
{
    totalTime++;
    if(funIndex==1)
        return five_uneven_peak_trap(x, dim);
    if(funIndex==2)
        return equal_maxima(x, dim);
    if(funIndex==3)
        return uneven_decreasing_maxima(x, dim);
    if(funIndex==4)
        return himmelblau(x, dim);
    if(funIndex==5)
        return six_hump_camel_back(x, dim);
    if(funIndex==6)
        return shubert(x, dim);
    if(funIndex==7)
        return vincent(x, dim);
    if(funIndex==8)
        return shubert(x, dim);
    if(funIndex==9)
        return vincent(x, dim);
    if(funIndex==10)
        return modified_rastrigin_all(x, dim);
    if(funIndex==11)
        return CF1(dim).evaluate(x);
    if(funIndex==12)
        return CF2(dim).evaluate(x);
    if(funIndex==13||funIndex==14||funIndex==16||funIndex==18)
        return CF3(dim).evaluate(x);
    if(funIndex==15||funIndex==17||funIndex==19||funIndex==20)
        return CF4(dim).evaluate(x);
    else
    {
        cout<<"wrong result,please check"<<endl;
        return -1000;
    }
}

//4. 计算两个体间的欧式距离
double caculateDis(int dim, double*a, double*b)
{
    double sum=0.0;
    for(int i=0;i<dim;++i)
        sum+=(a[i] - b[i])*(a[i] - b[i]);
    return sqrt(sum);
}

//5. a初始化边界

void initBound(int dim,double *UP,double*LOW)
{
    
    for (int i = 0; i < dim; ++i)
    {
        UP[i] = UB[i];
        LOW[i] = LB[i];
    }
    
    
}

//计算clpso相关参数
double calPc(int i,int numPop)
{
    double Pc = 0.0;
    Pc = 0.05 + 0.45*(1*exp(10*(i)*1.0/(numPop-1))-1)/(exp(10)-1);
    return Pc;
}


//6. 粒子类，包含定义以及更新方式
class Particle{
public:
    //构造函数：传递维度，小生境中心，Z值，上边界向量和下边界向量；同时初始化新产生粒子速度，并按相关属性计算适应值
    Particle(int dime, double*nicheC, double *Z, double* UP, double* LOW, double *pBestFit,double pBestCor[popNum][MAX_N],int Index )
    {
        dim = dime;
        for(int i=0;i<dime;++i)
            nicheCenter[i] = nicheC[i];
        for(int i=0;i<dime;++i){
            this->Z[i] = Z[i];
            cordinate[i]=rand_back2(Z[i]*(UP[i]-LOW[i]) + LOW[i], i);
            selfOptima[i] = cordinate[i];
        }
        for(int i=0;i<dime;++i)
            velocity[i] = 0;
        for(int i=0;i<dim;++i)
        {
            this->UP[i] = UP[i];
            this->LOW[i] = LOW[i];
        }
        for(int i=0;i<dim;++i)
        {
            if(cordinate[i]>UP[i])
                cordinate[i] = UP[i];
            if(cordinate[i]<LOW[i])
                cordinate[i] = LOW[i];
        }
        for(int i=0;i<dim;++i)
        {
            selfOptima[i] = cordinate[i];
        }
        for(int i=0;i<popNum;++i)
        {
            pbestFit[i] = pBestFit[i];
            for(int j=0;j<funDim;++j)
                pbestCor[i][j] = pBestCor[i][j];
        }
        fitness = funTest(dim, cordinate);
        selfOptimaFit = fitness;
        index = Index;
    }
    
    
    
    //按式(4)更新
    
    void explore(double*localBest,double *pBestFit,double pBestCor[popNum][MAX_N])
    {
        for(int i=0;i<popNum;++i)
        {
            pbestFit[i] = pBestFit[i];
            for(int j=0;j<funDim;++j)
                pbestCor[i][j] = pBestCor[i][j];
        }
        
        
        for(int i=0;i<funDim;++i)
        {
            
            velocity[i] = velocity[i]*w + Z1*r8_uniform_ab(0, 1)*(selfOptima[i] - cordinate[i]) +Z2*r8_uniform_ab(0, 1)*(localBest[i] - cordinate[i]);
            
        }
        
        for(int i=0;i<funDim;++i)
        {
            cordinate[i] += velocity[i];
        }
        for(int i=0;i<funDim;++i)
        {
            if(cordinate[i]>UP[i])
                cordinate[i] = UP[i];
            if(cordinate[i]<LOW[i])
                cordinate[i] = LOW[i];
        }
    }
    
    void exploreT(double*localBest,double *pBestFit,double pBestCor[popNum][MAX_N],int num)
    {
        for(int i=0;i<popNum;++i)
        {
            pbestFit[i] = pBestFit[i];
            for(int j=0;j<funDim;++j)
                pbestCor[i][j] = pBestCor[i][j];
        }
        
        
        for(int i=0;i<funDim;++i)
        {
            double Pc = calPc(index,num);
            double p = r8_uniform_ab(0, 1);
            if(p>Pc)
            {
                velocity[i] = velocity[i]*w +  Z1*r8_uniform_ab(0, 1)*(selfOptima[i] - cordinate[i]);
            }
            else
            {
                int rand1 = floor(r8_uniform_ab(0, 1)*num);
                int rand2 = floor(r8_uniform_ab(0, 1)*num);
                
                if(pbestFit[rand1]>pbestFit[rand2] && (pbestFit[rand1]-selfOptimaFit>0.0000000001||pbestFit[rand1]-selfOptimaFit<-0.0000000001))
                    velocity[i] = velocity[i]*w + Z1*r8_uniform_ab(0, 1)*(pbestCor[rand1][i] - cordinate[i]);
                else if(pbestFit[rand1]<=pbestFit[rand2] && (pbestFit[rand2]-selfOptimaFit>0.0000000001||pbestFit[rand2]-selfOptimaFit<-0.0000000001))
                    velocity[i] = velocity[i]*w + Z1*r8_uniform_ab(0, 1)*(pbestCor[rand2][i] - cordinate[i]);
                else
                {
                    while((pbestFit[rand2]-selfOptimaFit<0.0000000001&&pbestFit[rand2]-selfOptimaFit>-0.0000000001)&&rand2!=0)
                    {
                        rand2 = floor(r8_uniform_ab(0, 1)*num);
                    }
                    velocity[i] = velocity[i]*w + Z1*r8_uniform_ab(0, 1)*(pbestCor[rand2][i] - cordinate[i]);
                }
                
            }
        }
        
        
        for(int i=0;i<funDim;++i)
        {
            cordinate[i] += velocity[i];
        }
        for(int i=0;i<funDim;++i)
        {
            if(cordinate[i]>UP[i])
                cordinate[i] = UP[i];
            if(cordinate[i]<LOW[i])
                cordinate[i] = LOW[i];
        }
    }
    
    
    void modifyCor(int dim,double t)
    {
        cordinate[dim] = t;
        if(cordinate[dim]>UP[dim])
            cordinate[dim] = UP[dim];
        if(cordinate[dim]<LOW[dim])
            cordinate[dim] = LOW[dim];
        selfOptima[dim]=cordinate[dim];
    }
    
    void modifyfit()
    {
        fitness = funTest(dim, cordinate);
        selfOptimaFit = fitness;
    }
    
    
    //更新自身最优
    void updateParticle()
    {
        double temp = funTest(dim, cordinate);
        fitness = temp;
        if(temp>selfOptimaFit)
        {
            selfOptimaFit = temp;
            for(int i=0;i<dim;++i)
            {
                selfOptima[i]=cordinate[i];
            }
        }
        
    }
    
    double* getCordinate()
    {
        return cordinate;
    }
    
    double getSelfOptima()
    {
        return selfOptimaFit;
    }
    
    double *getSelfOptimaCor()
    {
        return selfOptima;
    }
    
    double getFitness()
    {
        return fitness;
    }
    
    int index;
    
    bool isModified = 0;
    
private:
    int dim;
    double w = 0.5;
    double Z1 = 2;
    double Z2 = 2;
    
    double Z[MAX_N];
    double UP[MAX_N];
    double LOW[MAX_N];
    double nicheCenter[MAX_N];
    
    double pbestFit[popNum];
    double pbestCor[popNum][MAX_N];
    double fitness;
    double selfOptimaFit;
    double selfOptima[MAX_N];      //pBest
    double cordinate[MAX_N];
    double velocity[MAX_N];
};


//7. 用于比较
class ComplessPop
{
public:
    bool operator()(Particle a, Particle b)
    {
        return a.getFitness() >b.getFitness();
    }
};

//8. 用于初始化小生境中心
double center[MAX_N];
double tZ[MAX_N];

double* initZ()
{
    for(int i=0;i<funDim;++i)
        tZ[i] = r8_uniform_ab(0, 1);
    return tZ;
}

//9. niche小生境类
class niche{
public:
    //构造函数: 用于传递维度，小生境大小，小生境半径，小生境中心，小生境Z值，小生境内部粒子更新次数;并计算小生境中心适应值
    niche(int dim, int size, double*center, double *Z, int time)
    {
        this->dim = dim;
        this->size = size;
        this->exploreTime = time;
        for(int i=0;i<dim;++i){
            this->Z[i] = Z[i];
            nicheCenter[i] = Z[i]*(UB[i]-LB[i])+LB[i];
            //          tempCenter[i] = center[i];
            UP[i] = UB[i];
            LOW[i] = LB[i];
        }
        fitness = funTest(dim, nicheCenter);
        tempFitness = -1000;
        for(int i=0;i<popNum;++i)
        {
            pbestFit[i] = -999999999;
            for(int j=0;j<dim;++j)
                pbestCor[i][j] = -9999;
        }
    }
    
    //小生境生成时，按大小初始化生成一定数量的粒子
    void createParticle()
    {
        for(int i=0;i<size;++i)
        {
            Particle temp(dim,nicheCenter,Z,UP,LOW,pbestFit,pbestCor,i);
            pop.push_back(temp);
        }
    }
    
    void sortPop()
    {
        sort(pop.begin(), pop.end(), ComplessPop());
    }
    
    void initialParticle()
    {
        sort(pop.begin(), pop.end(), ComplessPop());
        for(int i=pop.size()-1;i<pop.size();++i)
        {
            for(int k=0;k<pop.size();++k)
                pop[k].isModified = 0;
            for(int j=0;j<dim;++j)
            {
                pop[i].modifyCor(j, rand_back(nicheCenter[j], j));
                pop[i].isModified = 1;
            }
            
            pop[i].modifyfit();
        }
        
        
        nicheUpdate();
    }
    
    //按一定次数更新小生境内部粒子
    void explore()
    {
        for(int t=0;t<exploreTime;++t)
        {
            
            for(int i=0;i<pop.size();++i)
            {
                pop[i].explore(nicheCenter,pbestFit,pbestCor);
                pop[i].updateParticle();
            }
            nicheUpdate();
        }
    }
    
    /*
     void modifyCor(int dim,double cor)
     {
     nicheCenter[dim] = cor;
     }
     */
    
    void exploreT()
    {
        for(int t=0;t<exploreTime;++t)
        {
            if(calInnerDis()<disToSwitch)
            {
                for(int i=0;i<pop.size();++i)
                {
                    pop[i].exploreT(nicheCenter,pbestFit,pbestCor,pop.size());
                    pop[i].updateParticle();
                }
            }
            else
            {
                for(int i=0;i<pop.size();++i)
                {
                    pop[i].explore(nicheCenter,pbestFit,pbestCor);
                    pop[i].updateParticle();
                }
            }
            nicheUpdate();
        }
    }
    
    
    //更新小生境中心
    void nicheUpdate()
    {
        int index_MAX = 0;
        tempFitness = pop[0].getFitness();
        
        for(int i=0;i<pop.size();++i)
        {
            pbestFit[i] = pop[i].getSelfOptima();
            for(int j=0;j<dim;++j)
                pbestCor[i][j] = pop[i].getSelfOptimaCor()[j];
            
        }
        for(int i=1;i<size;++i)
        {
            if(pop[i].getFitness()>tempFitness)
            {
                index_MAX = i;
                tempFitness = pop[i].getFitness();
            }
        }
        if(tempFitness-fitness>0.00000000000001)
        {
            
            if(tempFitness-fitness>accuracy/10.0)
                suspendTime = 0;
            
            /*
             if(tempFitness-fitness>accuracy/100.0)
             suspendTime = 0;
             */
            fitness = tempFitness;
            for(int i=0;i<dim;++i)
            {
                nicheCenter[i] = pop[index_MAX].getCordinate()[i];
            }
        }
        else
            suspendTime++;
        
        
        
    }
    
    
    //将小生境种群压至规定大小，多退少补
    void allocatePOP(int newSize)
    {
        int excceedNum = (int)pop.size()-newSize;
        if(excceedNum<0)
        {
            for(int i=0;i<-excceedNum;++i)
            {
                Particle temp(dim,nicheCenter,Z,UP,LOW,pbestFit,pbestCor,i);
                pop.push_back(temp);
            }
        }
        
        else if(excceedNum>0)
        {
            sort(pop.begin(), pop.end(), ComplessPop());
            for(int i=0;i<pop.size();++i)
                pop[i].index = i;
            for(int i=0;i<excceedNum;++i)
            {
                pop.pop_back();
            }
        }
        size = (int)pop.size();
    }
    
    
    
    //更新Z值，并将该小生境全部map（仅小生境中心Z值起作用）
    void map(vector<niche>&b)
    {
        initBound(dim, UB, LB);
        double tempNicheCenter[MAX_N];
        for(int i=0;i<dim;++i)
        {
            tempNicheCenter[i] = Z[i]*(UB[i]-LB[i]) + LB[i];
            if(tempNicheCenter[i]>UP[i])
                tempNicheCenter[i] = UP[i];
            if(tempNicheCenter[i]<LOW[i])
                tempNicheCenter[i] = LOW[i];
        }
        for(int i=0;i<funDim;++i)
            circleMap(Z[i]);
        double newMap[MAX_N];
        for(int i=0;i<funDim;++i)
        {
            newMap[i] = r8_uniform_ab(0, 1);
        }
        niche temp(dim, 1,tempNicheCenter,newMap,exploreTime);
        temp.createParticle();
        temp.nicheUpdate();
        b.push_back(temp);
        
        
    }
    
    
    
    double calInnerDis()
    {
        double temp = 0;
        for(int i=0;i<pop.size();++i)
        {
            double temp2 = 0.0;
            for(int j=0;j<dim;++j)
                temp2 += (nicheCenter[j]-pop[i].getCordinate()[j])*(nicheCenter[j]-pop[i].getCordinate()[j]);
            temp+=sqrt(temp2);
        }
        temp /=pop.size();
        return temp;
    }
    
    
    double* getZ()
    {
        return Z;
    }
    
    double getFitness()
    {
        return fitness;
    }
    
    double* getCenter()
    {
        return nicheCenter;
    }
    
    int getSize()
    {
        return pop.size();
    }
    
    vector<Particle> getPop()
    {
        return pop;
    }
    
    int suspendTime = 0;
    bool isSuspended = 0;
    
private:
    int exploreTime;
    int dim;
    int size;
    
    double pbestCor[popNum][MAX_N];
    double pbestFit[popNum];
    double fitness;
    double tempFitness;
    double Z[MAX_N];
    double UP[MAX_N];
    double LOW[MAX_N];
    double nicheCenter[MAX_N];
    vector<Particle> pop;
};


class ComplessNiche
{
public:
    bool operator()(niche a, niche b)
    {
        return a.getFitness()>b.getFitness();
    }
};

class suspendSet
{
public:
    suspendSet (int a)
    {
        suspendFitness = INT_MIN;
        suspendNiche.clear();
    }
    bool addSuspendNiche(niche &a)
    {
        if(a.getFitness()-suspendFitness>0.0000000001&&a.suspendTime>=TimeToSuspendSet)
        {
            suspendFitness = a.getFitness();
            for(int i=suspendNiche.size()-1;i>=0;i--)
            {
                if(a.getFitness()-suspendNiche[i].getFitness()>accuracy/10)
                {
                    for(int k=i;k<suspendNiche.size()-1;++k)
                    {
                        suspendNiche[k] = suspendNiche[k+1];
                    }
                    suspendNiche.pop_back();
                }
            }
            a.isSuspended = true;
            suspendNiche.push_back(a);
            return false;
        }
        
        else if(a.getFitness()-suspendFitness>0.0000000001&&a.suspendTime<TimeToSuspendSet&&suspendNiche.size()!=0)
        {
            suspendFitness = a.getFitness();
            for(int i=suspendNiche.size()-1;i>=0;i--)
            {
                if(suspendNiche[i].getFitness()-a.getFitness()<-accuracy/10.0)
                {
                    
                    for(int k=i;k<suspendNiche.size()-1;++k)
                    {
                        suspendNiche[k] = suspendNiche[k+1];
                    }
                    suspendNiche.pop_back();
                    
                }
            }
            a.isSuspended = false;
            return false;
        }
        
        else if(a.getFitness() - suspendFitness>=-accuracy/10.0 && suspendNiche.size()!=0&&a.suspendTime>=TimeToSuspendSet)
        {
            suspendNiche.push_back(a);
            a.isSuspended = true;
        }
        
        else
            a.isSuspended = false;
        return true;
    }
    
    void setUpdate()
    {
        sort(suspendNiche.begin(),suspendNiche.end(),ComplessNiche());
        for(int i=0;i<suspendNiche.size();++i)
        {
            for(int j=i+1;j<suspendNiche.size();++j)
            {
                if(caculateDis(funDim, suspendNiche[i].getCenter(), suspendNiche[j].getCenter())<0.01)
                {
                    for(int k=j;k<suspendNiche.size()-1;++k)
                    {
                        suspendNiche[k] = suspendNiche[k+1];
                    }
                    suspendNiche.pop_back();
                    --j;
                }
            }
        }
    }
    
    double getFitness()
    {
        return suspendFitness;
    }
    vector<niche> getNiche()
    {
        return suspendNiche;
    }
    vector<niche> suspendNiche;
private:
    
    double suspendFitness = -10000;
};




double* initNicheCenter()
{
    for(int i=0;i<funDim;++i)
        center[i] = r8_uniform_ab(LB[i], UB[i]);
    return center;
}



double positive(double a)
{
    if(a>0)
        return a;
    else
        return -a;
}

double minDisBetNiche(vector<niche> Niche)
{
    double minDis = INT_MAX;
    for(int i=0;i<Niche.size();++i)
    {
        if(Niche[i].isSuspended==1)
            for(int j=i+1;j<Niche.size();++j)
            {
                if(Niche[i].isSuspended==1&&Niche[j].isSuspended==1)
                    if(caculateDis(funDim, Niche[i].getCenter(), Niche[j].getCenter())<minDis)
                        minDis =caculateDis(funDim, Niche[i].getCenter(), Niche[j].getCenter());
            }
    }
    return minDis;
}



int main(int argc, const char * argv[]) {
    double result = 0;
    double result2 = 0;
    double result3 = 0;
    double result4 = 0;
    double result5 = 0;
    initVal(funIndex);
    for(int dt=0;dt<30;++dt)
    {
        initVal(funIndex);
        setDisToSwitch();
        totalTime = 0;
        suspendSet SuspendSet(0);
        //初始化最初的initnumOfNiche个小生境,并对每个niche进行innerRunTime次内部更新
        srand((unsigned)time(0));
        vector<niche> Niche;                   //由所有niche组成的vector
        vector<double> Partion;                //以适应度值
        for(int i=0;i<initnumOfNiche;++i)
        {
            niche temp(funDim,initsizeOfNiche,initNicheCenter(),initZ(),innerRunTime);
            temp.createParticle();
            //temp.initialParticle();
            temp.exploreT();
            Niche.push_back(temp);
        }
        
        
        //对initnumOfNiche个小生境排序,为第一次division找到最差小生境
        sort(Niche.begin(), Niche.end(), ComplessNiche());
        
        
        for(int t=0;totalTime<maxFes;++t)
        {
            double partSum = 0;
            //将大小为popNum的种群划分到各小生境，以相应的适应度值为标准，Partion记录每个niche应当分配的粒子个数
            Partion.clear();
            int notSuspendNum = 0;
            for(int i=0;i<Niche.size();++i)
            {
                if(Niche[i].isSuspended == 0)
                {
                    partSum += Niche[i].getFitness() - Niche[Niche.size() - 1].getFitness()+0.0000000001;   //＝0.000000001是为了防止除0
                    notSuspendNum++;
                }
            }
            for(int i = 0;i<Niche.size();++i)
            {
                if(Niche[i].isSuspended == 0)
                {
                    if(ceil((Niche[i].getFitness() - Niche[Niche.size() - 1].getFitness() + 0.0000000001)/partSum*popNum)<10)
                        Partion.push_back(ceil((Niche[i].getFitness() - Niche[Niche.size() - 1].getFitness() + 0.0000000001)/partSum*popNum));
                    else
                        Partion.push_back(10);
                }
                else
                    Partion.push_back(0);
            }
            
            vector<niche> tempNiche ;    //用于把新生成的小生境暂时加入，之后会把其中新的小生境放回Niche
            for(int i=0;i<Niche.size();++i)
                tempNiche.push_back(Niche[i]);
            for(int i=0;i<Niche.size();++i)
            {
                if(Partion[i]==0)
                {
                    Niche[i].map(tempNiche);
                }
                else if(Partion[i]!=0)
                {
                    Niche[i].allocatePOP(Partion[i]);
                    Niche[i].exploreT();
                    Niche[i].map(tempNiche);
                    if(Niche[i].suspendTime>=TimeToSuspendSet||Niche[i].getFitness()-SuspendSet.getFitness()>-accuracy/10)
                    {
                        if(!SuspendSet.addSuspendNiche(Niche[i]))
                        {
                            for(int j=0;j<Niche.size();++j)
                            {
                                if(SuspendSet.getFitness()- Niche[j].getFitness()>accuracy/10&&j!=i)
                                    Niche[j].isSuspended = false;
                            }
                        }
                    }
                    Niche[i].initialParticle();
                }
            }
            
            SuspendSet.setUpdate();
            
            int exceed = tempNiche.size()-Niche.size();
            //把tempNiche中的新小生境放回niche
            while( exceed>0)
            {
                Niche.push_back(tempNiche[tempNiche.size()-1]);
                tempNiche.pop_back();
                exceed--;
            }
            
            
            
            //判断niche间是否过近，其中0.5那个是比较敏感的参数，在后期肯定是要把它函数化的
            double disOfDiagonal=0;
            for(int d=0;d<funDim;++d)
            {
                disOfDiagonal+=(UB[d]-LB[d])*(UB[d]-LB[d]);
            }
            disOfDiagonal = sqrt(disOfDiagonal);
            
            //同峰小生境删除
            sort(Niche.begin(), Niche.end(), ComplessNiche());
            for(int i=0;i<Niche.size();++i)
            {
                for(int j=i+1;j<Niche.size();++j)
                {
                    if(((caculateDis(funDim, Niche[i].getCenter(), Niche[j].getCenter())<disOfDiagonal*0.15))&&Niche[j].isSuspended==0)
                    {
                        double tempCor[MAX_N];
                        for(int d=0;d<funDim;++d)
                            tempCor[d] = (Niche[i].getCenter()[d]+Niche[j].getCenter()[d])/2.0;
                        double tempFit = funTest(funDim, tempCor);
                        if(!(Niche[i].getFitness()>tempFit&&Niche[j].getFitness()>tempFit))
                        {
                            
                            for(int k=j;k<Niche.size()-1;++k)
                            {
                                Niche[k] = Niche[k+1];
                            }
                            Niche.pop_back();
                            j--;
                        }
                    }
                    
                }
            }
            
            //较差小生境删除
            int activeNiche = 0;
            for(int i=0;i<Niche.size();++i)
                if(Niche[i].isSuspended == 0)
                    activeNiche++;
            //如果niche数量过大超高上限，则删除较差的多余niche
            if(activeNiche>MAX_NICHE_NUM)
            {
                int exceedNum = (int)activeNiche - MAX_NICHE_NUM;
                for(int i=0;i<exceedNum;++i)
                    Niche.pop_back();
            }
            
            
            //disToSwitch更新
            if(SuspendSet.getNiche().size()>1)
            {
                disToSwitch = min(minDisBetNiche(Niche)*1.0,disToSwitch);
            }
            /*
             for(int i=0;i<Niche.size();++i)
             {
             for(int j=0;j<funDim;++j)
             cout<<Niche[i].getCenter()[j]<<"   ";   //这里是坐标
             cout<<"result:"<<endl;
             cout<<Niche[i].getFitness()<<"  "<<Niche[i].getZ()[0]<<endl;
             }
             */
            
            
        }
        
        
        
        
        //用于测试
        for(int i=0;i<Niche.size();++i)
        {
            for(int j=i+1;j<Niche.size();++j)
            {
                if(caculateDis(funDim, Niche[i].getCenter(), Niche[j].getCenter())<rhO)
                {
                    for(int k=j;k<Niche.size()-1;++k)
                    {
                        Niche[k] = Niche[k+1];
                    }
                    Niche.pop_back();
                    --j;
                }
            }
        }
        //！！！！只需要修改这边就可以了，你先自己定义一个符合cec标准的关于解的结构体，然后把暂停集里面每个粒子和最后一次迭代的全部主粒子用你定义的结构体表示，作为api的参数吧
        int a=0;
        int b = 0;
        int c = 0;
        int d = 0;
        int e = 0;
        int f = 0;
        //这边是输出最后一次迭代的主粒子的有关信息
        
        for(int i=0;i<Niche.size();++i)
        {
            /*
             cout<<"result:"<<endl;
             
             if(Niche[i].isSuspended==0)
             cout<<Niche[i].calInnerDis()<<"    "<<Niche[i].getSize()<<"!!!"<<endl;
             */
            
            //cout<<Niche[i].getFitness()<<endl;
            
            if(maxFitness-Niche[i].getFitness()>0.00001)
                a++;
            if(maxFitness-Niche[i].getFitness()>0.0001)
                b++;
            if(maxFitness-Niche[i].getFitness()>0.001)
                c++;
            if(maxFitness-Niche[i].getFitness()>0.01)
                d++;
            if(maxFitness-Niche[i].getFitness()>0.1)
                f++;
            if(maxFitness-Niche[i].getFitness()<0.0001&&maxFitness-Niche[i].getFitness()>0.00001&&Niche[i].isSuspended==0)
                e++;
        }
        /*
         for(int i=0;i<Niche.size();++i)
         {
         for(int j=0;j<funDim;++j)
         cout<<Niche[i].getCenter()[j]<<"   ";   //这里是坐标
         cout<<"result:"<<endl;
         cout<<Niche[i].getFitness()<<"  "<<Niche[i].getZ()[0]<<endl;
         }
         */
        cout<<dt<<"th:"<<funDim<<endl;
        cout <<Niche.size()-a<<"  "<<Niche.size()-b<<"  "<<Niche.size()-c<<"  "<<Niche.size()-d<<"  "<<Niche.size()-f<<"  "<<e<<"  "<<totalTime<<endl;
        cout<<MAX_NICHE_NUM<<"  "<<minDisBetNiche(Niche)<<endl;
        
        //这里是输出暂停集的有关信息
        /*
         cout<<SuspendSet.getNiche().size()<<endl;
         for(int i=0;i<SuspendSet.getNiche().size();++i)
         {
         for(int j=0;j<funDim;++j)
         cout<<SuspendSet.getNiche()[i].getCenter()[j]<<"   ";  //这里是坐标
         cout<<"result:"<<endl;
         cout<<SuspendSet.getNiche()[i].getFitness()<<"  "<<SuspendSet.getNiche()[i].getZ()[0]<<endl;
         }*/
        for(int i=0;i<Niche.size();++i)
        {
            for(int j=0;j<funDim;++j)
                cout<<Niche[i].getCenter()[j]<<"   ";
            cout<<"result:"<<endl;
            cout<<Niche[i].getFitness()<<"  "<<Niche[i].getZ()[0]<<endl;
        }
        
        result+=Niche.size()-a;
        result2+=Niche.size()-b;
        result3+=Niche.size()-c;
        result4+=Niche.size()-d;
        result5+=Niche.size()-f;
    }
    cout<<result/30<<endl;
    cout<<result2/30<<endl;
    cout<<result3/30<<endl;
    cout<<result4/30<<endl;
    cout<<result5/30<<endl;
    return 0;
}
