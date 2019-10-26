#include "XAJModel.h"
XAJ::XAJ(void)//构造函数
{
	this->P = NULL;
	this->EM = NULL;
	this->W = NULL;
	this->R = NULL;
	this->RS = NULL;
	this->RSS = NULL;
	this->RG = NULL;
	this->QRS = NULL;
	this->QRSS = NULL;
	this->QRG = NULL;
	this->Q = NULL;
	this->WU = NULL;
	this->WL = NULL;
	this->WD = NULL;
	this->EU = NULL;
	this->EL = NULL;
	this->ED = NULL;
	this->E = NULL;
	this->RF = NULL;
}


XAJ::~XAJ(void)//析构函数
{
	delete[] this->P;
	delete[] this->EM;
	delete[] this->W;
	delete[] this->R;
	delete[] this->RS;
	delete[] this->RSS;
	delete[] this->RG;
	delete[] this->QRS;
	delete[] this->QRSS;
	delete[] this->QRG;
	delete[] this->Q;
	delete[] this->WU;
	delete[] this->WL;
	delete[] this->WD;
	delete[] this->EU;
	delete[] this->EL;
	delete[] this->ED;
	delete[] this->E;
	delete[] this->RF;
}

//初始化
void XAJ::Initialize(long steps, float Area, int deltaT, int IsRouting, char *inputfile, float WU0, float WL0, float WD0, float FR0, float S0, float QRSS0, float QRG0)
{
	/*
	参数说明
	steps 计算步长数
	Area 流域面积
	deltaT 计算步长（小时）
	IsRouting 是否汇流计算（1:是；0:否）
	inputfile 降雨蒸发输入文件
	WU0 上层初始含水量
	WL0 下层初始含水量
	WD0 深层初始含水量
	FR  流域初始自由水蓄水量小于等于SS的比值
	S 流域初始平均自由水蓄水容量
	QRSS0 流域出口壤中流初始流量 m3/s
	QRG0 流域出口地下径流初始流量 m3/s
	*/
	FILE * fp;
	this->steps = steps;
	//为模型参量分配内存，并初始化为0
	this->P = new float[this->steps]{ 0.0 };
	this->EM = new float[this->steps]{ 0.0 };
	this->W = new float[this->steps]{ 0.0 };
	this->R = new float[this->steps]{ 0.0 };
	this->RS = new float[this->steps]{ 0.0 };
	this->RSS = new float[this->steps]{ 0.0 };
	this->RG = new float[this->steps]{ 0.0 };
	this->QRS = new float[this->steps]{ 0.0 };
	this->QRSS = new float[this->steps]{ 0.0 };
	this->QRG = new float[this->steps]{ 0.0 };
	this->Q = new float[this->steps]{ 0.0 };
	this->WU = new float[this->steps]{ 0.0 };
	this->WL = new float[this->steps]{ 0.0 };
	this->WD = new float[this->steps]{ 0.0 };
	this->EU = new float[this->steps]{ 0.0 };
	this->EL = new float[this->steps]{ 0.0 };
	this->ED = new float[this->steps]{ 0.0 };
	this->E = new float[this->steps]{ 0.0 };
	this->RF = new float[this->steps]{ 0.0 };

	this->Area = Area;
	this->deltaT = deltaT;
	this->IsRouting = IsRouting;
	//打开输入文件，读入输入项（P和E）
	fp = fopen(inputfile, "r");
	for (int i = 0; i < this->steps; i++)
	{
		fscanf(fp, "%f %f", &(this->P[i]), &(this->EM[i]));
	}
	this->WU0 = WU0;
	this->WL0 = WL0;
	this->WD0 = WD0;
	this->FR0 = FR0;
	this->S0 = S0;
	this->QRSS0 = QRSS0;
	this->QRG0 = QRG0;
}
//参数设置（通过文件读入）
void XAJ::SetParameters(char *parfile)
{
	FILE * fp;
	fp = fopen(parfile, "r");
}

//参数设置（通过数组输入）
void XAJ::SetParameters(float *Parameters, float *UH)
{
	this->K = Parameters[0];
	this->C = Parameters[1];
	this->IMP = Parameters[2];
	this->WUM = Parameters[3];
	this->WLM = Parameters[4];
	this->WDM = Parameters[5];
	this->B = Parameters[6];
	this->SM = Parameters[7];
	this->EX = Parameters[8];
	this->KG = Parameters[9];
	this->KSS = Parameters[10];
	this->KKG = Parameters[11];
	this->KKSS = Parameters[12];
	this->n_unit = int(Parameters[13]);
	//部分参数可由输入参数计算得到
	this->WM = this->WUM + this->WLM + this->WDM;
	this->WMM = this->WM * (1.0 + this->B) / (1.0 - this->IMP);

	this->UH = new float[this->n_unit]{ 0.0 };
	for (int i = 0; i < this->n_unit; i++)
	{
		this->UH[i] = UH[i];
	}
}

//运行模型
void XAJ::RunModel(void)
{
	//首先定义步长内各个状态参量
	float PE_t;//净雨量
	float EP_t;//KC*EM
	float P_t;//降雨
	float R_t;//产流深
	float RB_t;//不透水面上产生的径流深度
	float RG_t;//地下径流深
	float RSS_t;//壤中流深度
	float RS_t;//地表径流深


	float S = S0; //自由水蓄水深
	float FS;//自由水蓄水量小于等于SS的面积
	float FR = FR0;//产流面积相对值，即FS与Area的比值
	float SSM; // (1+EX)*SM，流域上自由水蓄水量最大的某点的蓄量值
	float AU;//与自由蓄水量S对应的自由水蓄水容量曲线的纵坐标值
	float A;//土壤湿度为 W时土壤含水量折算成的径流深度（ mm）
	float E_t = 0.0;//蒸散发
	float EU_t = 0.0;//上层土壤蒸散发
	float EL_t = 0.0;//下层土壤蒸散发
	float ED_t = 0.0;//深层土壤蒸散发
	float WU_t = WU0;
	float WL_t = WL0;
	float WD_t = WD0;
	float W_t = WU_t + WL_t + WD_t;

	SSM = SM * (1 + EX);
	//FR0 = 1 - pow((1 - S0/ SSM), EX);

	int i; //代表步长下标
	for (i = 0; i < steps; i++)
	{
		/*蒸散发计算开始*/
		P_t = P[i] * (1 - IMP); //降在透水面的降雨量
		EP_t = K * EM[i];//降雨期间蒸发量
		if (P[i] > EP_t)
			RB_t = (P[i] - EP_t) * IMP;//RB是降在不透水面的降雨量
		else
			RB_t = 0.0;
		if ((WU_t + P_t) >= EP_t)//上层张力蓄水量足够
		{
			EU_t = EP_t;
			EL_t = 0;
			ED_t = 0;
		}
		else if ((WU_t + P_t) < EP_t) //上层张力蓄水量不够
		{
			EU_t = WU_t + P_t;//上层蒸发量为上层蓄水量+降雨量，上层变干
			EL_t = (EP_t - EU_t) * WL_t / WLM;
			if (EL_t < C * (EP_t - EU_t) && WL_t >= C * (EP_t - EU_t)) //
			{
				EL_t = C * (EP_t - EU_t);
				ED_t = 0;
			}
			else if (EL_t < C * (EP_t - EU_t) && WL_t < C * (EP_t - EU_t))//下层蓄量不够，触及深层
			{
				EL_t = WL_t;
				ED_t = C * (EP_t - EU_t) - EL_t;
			}
		}
		E_t = EU_t + EL_t + ED_t;
		PE_t = P_t - E_t;
		/*蒸散发计算结束*/

		/*子流域产流量计算开始*/
		if (PE_t <= 0) //净雨小于0，降雨全部蒸发
		{
			R_t = 0.0;//不产流
			W_t = W_t + PE_t;//更新含水量
		}
		else
		{

			A = WMM * (1 - pow((1.0 - W_t / WM), 1.0 / (1 + B)));
			// 土壤湿度折算净雨量 +降水后蒸发剩余雨量 <流域内最大含水容量
			if ((A + PE_t) < WMM)
			{
				R_t = PE_t + W_t + WM * pow((1 - (PE_t + A) / WMM), (1 + B)) - WM + RB_t;
			}
			// 土壤湿度折算净雨量 +降水后蒸发剩余雨量 >流域内最大含水容量
			else
			{
				// 流域内的产流深度计算
				R_t = PE_t + W_t - WM + RB_t;
			}
		}
		//三层蓄水量的计算：WU，WL，WD
		if (WU_t + P_t - EU_t - R_t <= WUM)//表层未达到蓄水容量
		{
			WU_t = WU_t + P_t - EU_t - R_t;
			WL_t = WL_t - EL_t;
			WD_t = WD_t - ED_t;
		}
		else
		{
			WU_t = WUM;//表层达到蓄水容量
			if (WL_t - EL_t + (WU_t + P_t - EU_t - R_t - WUM) < WLM)//下层未达到蓄水容量
			{
				WL_t = WL_t - EL_t + (WU_t + P_t - EU_t - R_t - WUM);
				WD_t = WD_t - ED_t;
			}
			else//下层达到蓄水容量
			{
				WL_t = WLM;
				if (WD_t - ED_t + WL_t - EL_t + (WU_t + P_t - EU_t - R_t - WUM) - WLM <= WDM)//深层未达到蓄水容量
					WD_t = WD_t - ED_t + WL_t - EL_t + (WU_t + P_t - EU_t - R_t - WUM) - WLM;
				else
					WD_t = WDM;
			}
		}
		/*子流域产流量计算结束*/

		/*三水源划分汇流计算*/
		if (PE_t > 0)//如果净雨大于0
		{
			FR = (R_t - RB_t) / PE_t;
			S = S0;
			AU = SSM * (1 - pow((1 - S / SM), 1 / (1 + EX)));
			if (PE_t + AU < SSM)
			{
				RS_t = FR * (PE_t + S - SM + SM * pow((1 - (PE_t + AU) / SSM), EX + 1));
				RSS_t = FR * KSS * (SM - SM * pow((1 - (PE_t + AU) / SSM), EX + 1));
				RG_t = FR * KG * (SM - SM * pow((1 - (PE_t + AU) / SSM), EX + 1));
				S0 = (1 - KSS - KG) * (SM - SM * pow((1 - (PE_t + AU) / SSM), EX + 1));
			}
			else if (PE_t + AU >= SSM)
			{
				RS_t = FR * (PE_t + S - SM);
				RSS_t = SM * KSS *FR;
				RG_t = SM * KG * FR;
				S0 = (1 - KSS - KG) * SM;
			}
			RS_t += RB_t;
			R_t = RS_t + RSS_t + RG_t;

			FR0 = FR;
		}
		else if (PE_t <= 0)
		{
			S = S0;
			FR = (1 - pow((1 - W_t / WM), B / (1 + B)));
			//RSS_t = 0.0;
			//RG_t = 0.0;
			RSS_t = S * KSS * FR;
			RG_t = S * KG * FR;
			RS_t = RB_t;
			R_t = RS_t + RSS_t + RG_t;

			S0 = S * (1 - KSS - KG);
			FR0 = FR;
		}
		/*三水源划分汇流计算结束*/
		//状态量保存
		this->E[i] = E_t;
		this->EU[i] = EU_t;
		this->EL[i] = EL_t;
		this->ED[i] = ED_t;
		this->W[i] = W_t;
		this->WU[i] = WU_t;
		this->WL[i] = WL_t;
		this->WD[i] = WD_t;
		this->RG[i] = RG_t;
		this->RS[i] = RS_t;
		this->RSS[i] = RSS_t;
		this->R[i] = R_t;
	}
}
//保存结果
void  XAJ::SaveOutput(char *output)
{
	int i;
	FILE *fp;
	if ((fp = fopen(output, "w")) == NULL)
	{
		printf("Can not create output file!\n");
		return;
	}
	fprintf(fp, "Timestep, E(mm), EU(mm), EL(mm), ED(mm), W(mm), WU(mm), WL(mm), WD(mm), R(mm), RS(mm), RSS(mm), RG(mm), Q(m3/d), QS(m3/d), QSS(m3/d), QG(m3/d)\n");
	for (i = 0; i < this->steps; i++)
	{
		fprintf(fp, "%d, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf, %9.3lf\n",
			i, this->E[i],this->EU[i], this->EL[i], this->ED[i],
			this->W[i], this->WU[i], this->WL[i], this->WD[i],
			this->R[i], this->RS[i], this->RSS[i], this->RG[i],
			this->Q[i], this->QRS[i], this->QRSS[i], this->QRG[i]);
	}
	// 当前步长的总产流径流深度fclose(fp);
}
//计算产流
void  XAJ::Routing()
{
	/*地面径流汇流计算：单位线法*/
	int i, j;
	float B[10000] = { 0.0 }; //开个大数组用于存储结果
	for (i = 0; i < this->steps; i++)
	{
		for (j = 0; j < this->n_unit; j++)
		{
			B[i + j] += this->RS[i] * this->UH[j];
		}
	}
	for (i = 0; i < steps; i++)
	{
		this->QRS[i] = B[i];
		// printf("QRS %f\n", this->QRS[i]);
	}
	/*壤中流汇流计算 :线性水库*/
	this->QRSS[0] = this->QRSS0;
	for (i = 1; i < this->steps; i++)
	{
		this->QRSS[i] = this->KKSS * this->QRSS[i - 1] + (1 - this->KKSS) * this->RSS[i] * this->Area / (3.6 * this->deltaT);
		// printf("QRSS %f\n", this->QRSS[i]);
	}
	/* 地下径流汇流计算 :线性水库*/
	this->QRG[0] = this->QRG0;
	for (int i = 1; i < this->steps; i++)
	{
		this->QRG[i] = this->KKG * this->QRG[i - 1] + (1 - this->KKG) * this->RG[i] * this->Area / (3.6 * this->deltaT);
		// printf("QRG %f\n", this->QRG[i]);
	}
	/*单元面积总入流计算*/
	for (i = 0; i < this->steps; i++)
	{
		this->Q[i] = this->QRS[i] + this->QRSS[i] + this->QRG[i];
		// printf("i:%d, Q: %.3f\n", i, this->Q[i]);
	}

}