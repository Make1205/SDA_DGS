#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <bits/stdc++.h>
#include <NTL/RR.h>
#include <random>

using namespace std;
using namespace NTL;

NTL_CLIENT
Vec<ZZ> SDAT;
Vec<RR> prob,probSDA,probfull,probtrun;
Vec<ZZ> probS;
Vec<RR> probSD;


void Frodo640Test()
{
    int s=12;

    int n=10000000;
    Vec<ZZ> OriT,SDATable;
    ifstream file1("../FrodoDist/Frodo640FP.txt");
    string str_num;
    for(int i=0;i<=s;i++)
    {
        getline(file1, str_num); // 读取文件中的字符串
        OriT.append(conv<ZZ>(str_num.c_str()));
    }
    for(int i=1;i<=s;i++)
    {
        OriT[i]+=OriT[i-1];
    }
    ifstream file2("../FrodoDist/Frodo640SDA.txt");
    for(int i=0;i<=s;i++)
    {
        getline(file2, str_num); // 读取文件中的字符串
        SDATable.append(conv<ZZ>(str_num.c_str()));
        // probSDA.append(to_RR(conv<ZZ>(str_num.c_str())));
    }
    for(int i=1;i<=s;i++)
    {
        SDATable[i]+=SDATable[i-1];
    }
    
    ZZ sgn;
    ZZ q=OriT[s]+ZZ(1);
    ZZ rnd1,rnd2;
    double t1,t2;
    int res=0;
    clock_t start_time1 = clock();
        for(int i=0;i<n;i++)
        {
            RandomBits(sgn,1);
            rnd1=RandomBnd(q);
            res=0;
            for(int i=0;i<=s;i++)
            {
                res += (rnd1 > OriT[i]);
            }
            if(sgn==ZZ(0)) res=-res; 
        }
    clock_t end_time1 = clock();
    t1 = ((double)(end_time1 - start_time1) / CLOCKS_PER_SEC );
    
    q=SDATable[s];
    // cout<<q<<endl;

    clock_t start_time2 = clock();
    for(int i=0;i<n;i++)
    {
        RandomBits(sgn,1);
        rnd2 = RandomBnd(q);
        res=0;
        for(int i=0;i<=s-1;i++)
        {
            res += (rnd2 > SDATable[i]);
        }
        if(sgn==ZZ(0)) res=-res; 
    }
    clock_t end_time2 = clock();
    t2 = ((double)(end_time2 - start_time2) / CLOCKS_PER_SEC );

    long mem0=0;
    for(int i=0;i<=s;i++)
    {
        mem0 += NumBits(OriT[i])* ( i >= 0 & OriT[i+1]!=OriT[i]);
    }
    long memS=0;
    for(int i=0;i<=s;i++)
    {
        memS += NumBits(SDATable[i])* ( i >= 0 & SDATable[i+1]!=SDATable[i]);
    }

    cout<<"Frodo640:"<<endl;
    cout<<"Speed:"<<endl;
    cout<<"Frodo640_Original_CDT: "<<to_RR(n)/t1<<" Samples/s"<<endl;
    cout<<"Frodo640_  SDA  _ CDT: "<<to_RR(n)/t2<<" Samples/s"<<endl;
    cout<<"Memory:"<<endl;
    cout<<"Frodo640_Original_CDT: "<<mem0<<" bits"<<endl;
    cout<<"Frodo640_  SDA  _ CDT: "<<memS<<" bits"<<endl;
    cout<<endl;
}


void Frodo976Test()
{

    int s=10;
    int n=10000000;
    Vec<ZZ> OriT,SDATable;
    ifstream file1("../FrodoDist/Frodo976FP.txt");
    string str_num;
    for(int i=0;i<=s;i++)
    {
        getline(file1, str_num); // 读取文件中的字符串
        OriT.append(conv<ZZ>(str_num.c_str()));
    }
    for(int i=1;i<=s;i++)
    {
        OriT[i]+=OriT[i-1];
    }

    ifstream file2("../FrodoDist/Frodo976SDA.txt");
    for(int i=0;i<=s;i++)
    {
        getline(file2, str_num); // 读取文件中的字符串
        SDATable.append(conv<ZZ>(str_num.c_str()));
        // probSDA.append(to_RR(conv<ZZ>(str_num.c_str())));
    }
    for(int i=1;i<=s;i++)
    {
        SDATable[i]+=SDATable[i-1];
    }

    ZZ sgn;
    ZZ q=OriT[s]+ZZ(1);
    ZZ rnd1,rnd2;
    double t1,t2;
    int res=0;
    clock_t start_time1 = clock();
        for(int i=0;i<n;i++)
        {
            RandomBits(sgn,1);
            rnd1=RandomBnd(q);
            res=0;
            for(int i=0;i<=s;i++)
            {
                res += (rnd1 > OriT[i]);
            }
            if(sgn==ZZ(0)) res=-res; 
        }
    clock_t end_time1 = clock();
    t1 = ((double)(end_time1 - start_time1) / CLOCKS_PER_SEC );
    q=SDATable[s];
    clock_t start_time2 = clock();
    for(int i=0;i<n;i++)
    {
        RandomBits(sgn,1);
        rnd2 = RandomBnd(q);
        res=0;
        for(int i=0;i<=s-1;i++)
        {
            res += (rnd2 > SDATable[i]);
        }
        if(sgn==ZZ(0)) res=-res; 
    }
    clock_t end_time2 = clock();
    t2 = ((double)(end_time2 - start_time2) / CLOCKS_PER_SEC );

    long mem0=0;
    for(int i=0;i<=s;i++)
    {
        mem0 += NumBits(OriT[i]);
    }
    long memS=0;
    for(int i=0;i<=s-1;i++)
    {
        memS += NumBits(SDATable[i]);
    }

    cout<<"Frodo976:"<<endl;
    cout<<"Speed:"<<endl;
    cout<<"Frodo976_Original_CDT: "<<to_RR(n)/t1<<" Samples/s"<<endl;
    cout<<"Frodo976_  SDA  _ CDT: "<<to_RR(n)/t2<<" Samples/s"<<endl;
    cout<<"Memory:"<<endl;
    cout<<"Frodo976_Original_CDT: "<<mem0<<" bits"<<endl;
    cout<<"Frodo976_  SDA  _ CDT: "<<memS<<" bits"<<endl;
    cout<<endl;
}

void Frodo1344Test()
{
    int s=6;
    int n=10000000;
    Vec<ZZ> OriT,SDATable;
    ifstream file1("../FrodoDist/Frodo1344FP.txt");
    string str_num;
    for(int i=0;i<=s;i++)
    {
        getline(file1, str_num); // 读取文件中的字符串
        OriT.append(conv<ZZ>(str_num.c_str()));
    }
    for(int i=1;i<=s;i++)
    {
        OriT[i]+=OriT[i-1];
    }
    ifstream file2("../FrodoDist/Frodo1344SDA.txt");
    for(int i=0;i<=s;i++)
    {
        getline(file2, str_num); // 读取文件中的字符串
        SDATable.append(conv<ZZ>(str_num.c_str()));
        // probSDA.append(to_RR(conv<ZZ>(str_num.c_str())));
    }
    for(int i=1;i<=s;i++)
    {
        SDATable[i]+=SDATable[i-1];
    }

    ZZ sgn;
    ZZ q=OriT[s]+ZZ(1);
    // cout<<q<<endl;
    ZZ rnd1,rnd2;
    double t1,t2;
    int res=0;
    // cout<<q<<endl;

    clock_t start_time1 = clock();
        for(int i=0;i<n;i++)
        {
            RandomBits(sgn,1);
            rnd1=RandomBnd(q);
            res=0;
            for(int i=0;i<=s;i++)
            {
                res += (rnd1 > OriT[i]);
            }
            if(sgn==ZZ(0)) res=-res; 
        }
    clock_t end_time1 = clock();
    t1 = ((double)(end_time1 - start_time1) / CLOCKS_PER_SEC );
    
    q=SDATable[s];
    // cout<<q<<endl;

    clock_t start_time2 = clock();
    for(int i=0;i<n;i++)
    {
        RandomBits(sgn,1);
        rnd2 = RandomBnd(q);
        res=0;
        for(int i=0;i<=s;i++)
        {
            res += (rnd2 > SDATable[i]);
        }
        if(sgn==ZZ(0)) res=-res; 
    }
    clock_t end_time2 = clock();
    t2 = ((double)(end_time2 - start_time2) / CLOCKS_PER_SEC );


    long mem0=0;
    for(int i=0;i<=s;i++)
    {
        mem0 += NumBits(OriT[i]);
    }
    long memS=0;
    for(int i=0;i<=s;i++)
    {
        memS += NumBits(SDATable[i]);
    }
    cout<<"Frodo1344:"<<endl;
    cout<<"Speed:"<<endl;
    cout<<"Frodo1344_Original_CDT: "<<to_RR(n)/t1<<" Samples/s"<<endl;
    cout<<"Frodo1344_  SDA  _ CDT: "<<to_RR(n)/t2<<" Samples/s"<<endl;
    cout<<"Memory:"<<endl;
    cout<<"Frodo1344_Original_CDT: "<<mem0<<" bits"<<endl;
    cout<<"Frodo1344_  SDA  _ CDT: "<<memS<<" bits"<<endl;
    cout<<endl;
}


void FalconTest()
{
    int s=18;
    int n=10000000;
    Vec<ZZ> OriT,SDATable;
    ifstream file1("../FalconDist/numberFP.txt");
    string str_num;
    for(int i=0;i<=s;i++)
    {
        getline(file1, str_num); // 读取文件中的字符串
        OriT.append(conv<ZZ>(str_num.c_str()));
    }
    for(int i=1;i<=s;i++)
    {
        OriT[i]+=OriT[i-1];
    }
    ifstream file2("../FalconDist/numberSDA.txt");
    for(int i=0;i<=s;i++)
    {
        getline(file2, str_num); // 读取文件中的字符串
        SDATable.append(conv<ZZ>(str_num.c_str()));
        // probSDA.append(to_RR(conv<ZZ>(str_num.c_str())));
    }
    for(int i=1;i<=s;i++)
    {
        SDATable[i]+=SDATable[i-1];
    }

    ZZ sgn;
    ZZ q=OriT[s];
    // cout<<q<<endl;
    ZZ rnd1,rnd2;
    double t1,t2;
    int res=0;
    // cout<<q<<endl;

    clock_t start_time1 = clock();
        for(int i=0;i<n;i++)
        {
            RandomBits(sgn,1);
            rnd1=RandomBnd(q);
            res=0;
            for(int i=0;i<=s;i++)
            {
                res += (rnd1 > OriT[i]);
            }
            if(sgn==ZZ(0)) res=-res; 
        }
    clock_t end_time1 = clock();
    t1 = ((double)(end_time1 - start_time1) / CLOCKS_PER_SEC );
    
    q=SDATable[s];
    // cout<<q<<endl;

    clock_t start_time2 = clock();
    for(int i=0;i<n;i++)
    {
        RandomBits(sgn,1);
        rnd2 = RandomBnd(q);
        res=0;
        for(int i=0;i<=s;i++)
        {
            res += (rnd2 > SDATable[i]);
        }
        if(sgn==ZZ(0)) res=-res; 
    }
    clock_t end_time2 = clock();
    t2 = ((double)(end_time2 - start_time2) / CLOCKS_PER_SEC );


    long mem0=0;
    for(int i=0;i<=s;i++)
    {
        mem0 += NumBits(OriT[i]);
    }
    long memS=0;
    for(int i=0;i<=s;i++)
    {
        memS += NumBits(SDATable[i]);
    }
    cout<<"Falcon:"<<endl;
    cout<<"Speed:"<<endl;
    cout<<"Falcon_Original_CDT: "<<to_RR(n)/t1<<" Samples/s"<<endl;
    cout<<"Falcon_  SDA  _ CDT: "<<to_RR(n)/t2<<" Samples/s"<<endl;
    cout<<"Memory:"<<endl;
    cout<<"Falcon_Original_CDT: "<<mem0<<" bits"<<endl;
    cout<<"Falcon_  SDA  _ CDT: "<<memS<<" bits"<<endl;
    cout<<endl;
}


int main()
{
    // cout<<"Speed test:"<<endl;
    Frodo640Test();
    Frodo976Test();
    Frodo1344Test();
    FalconTest();
    cout<<""<<endl;
    return 0;
}