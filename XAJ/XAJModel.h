#pragma once
#include "iostream"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
//�°���ģ��
class XAJ
{
private://ģ�Ͳ���
		//��һ�ࣺ��ɢ���������
	float K;//������ɢ��������������������֮��
	float WUM; //����ƽ���ϲ���ˮ��������λmm
	float WLM; //����ƽ���²���ˮ��������λmm
	float WDM; //����ƽ�������ˮ��������λmm
	float C;   //�����ɢ��ϵ��

			   //�ڶ��ࣺ�����������
	float IMP;//��͸ˮ�������
	float WM; //����ƽ����ˮ������ָ����ˮ��(WM=WUM+WLM+WDM��
	float B; //��ˮ��������ָ��

			 //�����ࣺˮԴ���ֲ���
	float SM; //������ˮ�������������������mm��
	float EX;//����ˮ��ˮ����~����ֲ�����ָ��
	float KSS; //�������ճ���ϵ��
	float KG;//����ˮ�ճ���ϵ��

			 //�����ࣺ���������еĲ���
	float KKSS; //����������ϵ��
	float KKG;//����ˮ����ˮϵ��
	float *UH; //����εĵر�����λ���ݱ�
	float KE;  //��Ԫ�Ӷε���˹����Kֵ
	float XE;  //��Ԫ�Ӷε���˹����Xֵ

			   //��������
	float WMM; //�����������ˮ��������λmm
	float Area;//�����������λkm2
	int deltaT;//����Сʱ��
	int IsRouting;//�Ƿ���кӵ���������
				  //һЩ��ʼ״̬����
	int n_unit;//����ε�λ��ʱ����
	float WU0;
	float WL0;
	float WD0;
	float FR0;
	float S0;
	float QRSS0;
	float QRG0;

public://ģ�Ͳ���
	   //������
	float *P;//����ˮ���λmm
	float *EM;//�����������λmm
	long steps;//���㲽����
			   //�����
	float *W;//����ʪ��
	float *R;//���������λmm
	float *RS;//�ر������λmm
	float *RSS;//���������λmm
	float *RG;//���¾������λmm
	float *QRS;//������ڵر�����
	float *QRSS;//���������������
	float *QRG;//������ڵ��¾�����
	float *Q;//�������������
	float U;//

			//������ز���
	float *WU;//�ϲ�����ʪ��
	float *WL;//�²�����ʪ��
	float *WD;//�������ʪ��
			  //������ز���
	float *EU;//�ϲ�����������(mm)
	float *EL;//�²�����������(mm)
	float *ED;//�������������(mm)
	float *E;//��������(mm)
			 //������ز���
	float *RF;

	XAJ(void);//���캯��
	~XAJ(void);//��������
			   //ģ�ͳ�ʼ��
	void Initialize(long steps, float Area, int deltaT, int IsRouting, char *inputfile, float WU0, float WL0, float WD0, float FR0, float S0, float QRSS0, float QRG0);
	//����ģ�Ͳ���
	void SetParameters(char *parfile);
	void SetParameters(float *Parameters, float *UH);//����
													 //����ģ��
	void RunModel(void);
	//������
	void SaveOutput(char *output);
	//�������
	void Routing();
};
