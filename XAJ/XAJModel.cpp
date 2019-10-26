#include "XAJModel.h"
XAJ::XAJ(void)//���캯��
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


XAJ::~XAJ(void)//��������
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

//��ʼ��
void XAJ::Initialize(long steps, float Area, int deltaT, int IsRouting, char *inputfile, float WU0, float WL0, float WD0, float FR0, float S0, float QRSS0, float QRG0)
{
	/*
	����˵��
	steps ���㲽����
	Area �������
	deltaT ���㲽����Сʱ��
	IsRouting �Ƿ�������㣨1:�ǣ�0:��
	inputfile �������������ļ�
	WU0 �ϲ��ʼ��ˮ��
	WL0 �²��ʼ��ˮ��
	WD0 ����ʼ��ˮ��
	FR  �����ʼ����ˮ��ˮ��С�ڵ���SS�ı�ֵ
	S �����ʼƽ������ˮ��ˮ����
	QRSS0 ���������������ʼ���� m3/s
	QRG0 ������ڵ��¾�����ʼ���� m3/s
	*/
	FILE * fp;
	this->steps = steps;
	//Ϊģ�Ͳ��������ڴ棬����ʼ��Ϊ0
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
	//�������ļ������������P��E��
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
//�������ã�ͨ���ļ����룩
void XAJ::SetParameters(char *parfile)
{
	FILE * fp;
	fp = fopen(parfile, "r");
}

//�������ã�ͨ���������룩
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
	//���ֲ������������������õ�
	this->WM = this->WUM + this->WLM + this->WDM;
	this->WMM = this->WM * (1.0 + this->B) / (1.0 - this->IMP);

	this->UH = new float[this->n_unit]{ 0.0 };
	for (int i = 0; i < this->n_unit; i++)
	{
		this->UH[i] = UH[i];
	}
}

//����ģ��
void XAJ::RunModel(void)
{
	//���ȶ��岽���ڸ���״̬����
	float PE_t;//������
	float EP_t;//KC*EM
	float P_t;//����
	float R_t;//������
	float RB_t;//��͸ˮ���ϲ����ľ������
	float RG_t;//���¾�����
	float RSS_t;//���������
	float RS_t;//�ر�����


	float S = S0; //����ˮ��ˮ��
	float FS;//����ˮ��ˮ��С�ڵ���SS�����
	float FR = FR0;//����������ֵ����FS��Area�ı�ֵ
	float SSM; // (1+EX)*SM������������ˮ��ˮ������ĳ�������ֵ
	float AU;//��������ˮ��S��Ӧ������ˮ��ˮ�������ߵ�������ֵ
	float A;//����ʪ��Ϊ Wʱ������ˮ������ɵľ�����ȣ� mm��
	float E_t = 0.0;//��ɢ��
	float EU_t = 0.0;//�ϲ�������ɢ��
	float EL_t = 0.0;//�²�������ɢ��
	float ED_t = 0.0;//���������ɢ��
	float WU_t = WU0;
	float WL_t = WL0;
	float WD_t = WD0;
	float W_t = WU_t + WL_t + WD_t;

	SSM = SM * (1 + EX);
	//FR0 = 1 - pow((1 - S0/ SSM), EX);

	int i; //�������±�
	for (i = 0; i < steps; i++)
	{
		/*��ɢ�����㿪ʼ*/
		P_t = P[i] * (1 - IMP); //����͸ˮ��Ľ�����
		EP_t = K * EM[i];//�����ڼ�������
		if (P[i] > EP_t)
			RB_t = (P[i] - EP_t) * IMP;//RB�ǽ��ڲ�͸ˮ��Ľ�����
		else
			RB_t = 0.0;
		if ((WU_t + P_t) >= EP_t)//�ϲ�������ˮ���㹻
		{
			EU_t = EP_t;
			EL_t = 0;
			ED_t = 0;
		}
		else if ((WU_t + P_t) < EP_t) //�ϲ�������ˮ������
		{
			EU_t = WU_t + P_t;//�ϲ�������Ϊ�ϲ���ˮ��+���������ϲ���
			EL_t = (EP_t - EU_t) * WL_t / WLM;
			if (EL_t < C * (EP_t - EU_t) && WL_t >= C * (EP_t - EU_t)) //
			{
				EL_t = C * (EP_t - EU_t);
				ED_t = 0;
			}
			else if (EL_t < C * (EP_t - EU_t) && WL_t < C * (EP_t - EU_t))//�²������������������
			{
				EL_t = WL_t;
				ED_t = C * (EP_t - EU_t) - EL_t;
			}
		}
		E_t = EU_t + EL_t + ED_t;
		PE_t = P_t - E_t;
		/*��ɢ���������*/

		/*��������������㿪ʼ*/
		if (PE_t <= 0) //����С��0������ȫ������
		{
			R_t = 0.0;//������
			W_t = W_t + PE_t;//���º�ˮ��
		}
		else
		{

			A = WMM * (1 - pow((1.0 - W_t / WM), 1.0 / (1 + B)));
			// ����ʪ�����㾻���� +��ˮ������ʣ������ <���������ˮ����
			if ((A + PE_t) < WMM)
			{
				R_t = PE_t + W_t + WM * pow((1 - (PE_t + A) / WMM), (1 + B)) - WM + RB_t;
			}
			// ����ʪ�����㾻���� +��ˮ������ʣ������ >���������ˮ����
			else
			{
				// �����ڵĲ�����ȼ���
				R_t = PE_t + W_t - WM + RB_t;
			}
		}
		//������ˮ���ļ��㣺WU��WL��WD
		if (WU_t + P_t - EU_t - R_t <= WUM)//���δ�ﵽ��ˮ����
		{
			WU_t = WU_t + P_t - EU_t - R_t;
			WL_t = WL_t - EL_t;
			WD_t = WD_t - ED_t;
		}
		else
		{
			WU_t = WUM;//���ﵽ��ˮ����
			if (WL_t - EL_t + (WU_t + P_t - EU_t - R_t - WUM) < WLM)//�²�δ�ﵽ��ˮ����
			{
				WL_t = WL_t - EL_t + (WU_t + P_t - EU_t - R_t - WUM);
				WD_t = WD_t - ED_t;
			}
			else//�²�ﵽ��ˮ����
			{
				WL_t = WLM;
				if (WD_t - ED_t + WL_t - EL_t + (WU_t + P_t - EU_t - R_t - WUM) - WLM <= WDM)//���δ�ﵽ��ˮ����
					WD_t = WD_t - ED_t + WL_t - EL_t + (WU_t + P_t - EU_t - R_t - WUM) - WLM;
				else
					WD_t = WDM;
			}
		}
		/*������������������*/

		/*��ˮԴ���ֻ�������*/
		if (PE_t > 0)//����������0
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
		/*��ˮԴ���ֻ����������*/
		//״̬������
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
//������
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
	// ��ǰ�������ܲ����������fclose(fp);
}
//�������
void  XAJ::Routing()
{
	/*���澶���������㣺��λ�߷�*/
	int i, j;
	float B[10000] = { 0.0 }; //�������������ڴ洢���
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
	/*�������������� :����ˮ��*/
	this->QRSS[0] = this->QRSS0;
	for (i = 1; i < this->steps; i++)
	{
		this->QRSS[i] = this->KKSS * this->QRSS[i - 1] + (1 - this->KKSS) * this->RSS[i] * this->Area / (3.6 * this->deltaT);
		// printf("QRSS %f\n", this->QRSS[i]);
	}
	/* ���¾����������� :����ˮ��*/
	this->QRG[0] = this->QRG0;
	for (int i = 1; i < this->steps; i++)
	{
		this->QRG[i] = this->KKG * this->QRG[i - 1] + (1 - this->KKG) * this->RG[i] * this->Area / (3.6 * this->deltaT);
		// printf("QRG %f\n", this->QRG[i]);
	}
	/*��Ԫ�������������*/
	for (i = 0; i < this->steps; i++)
	{
		this->Q[i] = this->QRS[i] + this->QRSS[i] + this->QRG[i];
		// printf("i:%d, Q: %.3f\n", i, this->Q[i]);
	}

}