#pragma once
#include "iostream"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
//新安江模型
class XAJ
{
private://模型参数
		//第一类：蒸散发计算参数
	float K;//流域蒸散发能力与蒸发皿蒸发量之比
	float WUM; //流域平均上层蓄水容量，单位mm
	float WLM; //流域平均下层蓄水容量，单位mm
	float WDM; //流域平均深层蓄水容量，单位mm
	float C;   //深层蒸散发系数

			   //第二类：产流计算参数
	float IMP;//不透水面积比重
	float WM; //流域平均蓄水容量（指张力水）(WM=WUM+WLM+WDM）
	float B; //蓄水容量曲线指数

			 //第三类：水源划分参数
	float SM; //自由蓄水库容量（即最大蓄量，mm）
	float EX;//自由水蓄水容量~面积分布曲线指数
	float KSS; //壤中流日出流系数
	float KG;//地下水日出流系数

			 //第四类：汇流计算中的参数
	float KKSS; //壤中流消退系数
	float KKG;//地下水日退水系数
	float *UH; //无因次的地表径流单位线纵表
	float KE;  //单元河段的马斯京根K值
	float XE;  //单元河段的马斯京根X值

			   //其他参数
	float WMM; //流域内最大蓄水容量，单位mm
	float Area;//流域面积，单位km2
	int deltaT;//步长小时数
	int IsRouting;//是否进行河道汇流计算
				  //一些初始状态参量
	int n_unit;//无因次单位线时段数
	float WU0;
	float WL0;
	float WD0;
	float FR0;
	float S0;
	float QRSS0;
	float QRG0;

public://模型参量
	   //输入项
	float *P;//流域降水深，单位mm
	float *EM;//流域蒸发深，单位mm
	long steps;//计算步长数
			   //输出项
	float *W;//土壤湿度
	float *R;//流域径流深，单位mm
	float *RS;//地表径流深，单位mm
	float *RSS;//壤中流深，单位mm
	float *RG;//地下径流深，单位mm
	float *QRS;//流域出口地表径流量
	float *QRSS;//流域出口壤中流量
	float *QRG;//流域出口地下径流量
	float *Q;//流域出口总流量
	float U;//

			//土壤相关参量
	float *WU;//上层土壤湿度
	float *WL;//下层土壤湿度
	float *WD;//深层土壤湿度
			  //蒸发相关参量
	float *EU;//上层土壤蒸发量(mm)
	float *EL;//下层土壤蒸发量(mm)
	float *ED;//深层土壤蒸发量(mm)
	float *E;//总蒸发量(mm)
			 //径流相关参量
	float *RF;

	XAJ(void);//构造函数
	~XAJ(void);//析构函数
			   //模型初始化
	void Initialize(long steps, float Area, int deltaT, int IsRouting, char *inputfile, float WU0, float WL0, float WD0, float FR0, float S0, float QRSS0, float QRG0);
	//设置模型参数
	void SetParameters(char *parfile);
	void SetParameters(float *Parameters, float *UH);//重载
													 //运行模型
	void RunModel(void);
	//保存结果
	void SaveOutput(char *output);
	//计算产流
	void Routing();
};
