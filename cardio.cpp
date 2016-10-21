#include <armadillo>
#include <iomanip>
#include <sstream>
#include <iostream>
#include "wavelet2d.h"
#include "fftw3.h"
#include "cardio.h"
#include "Storm.h"

using namespace std;
using namespace arma;
using namespace storm;


/**************************************************************
*函数名:		diff
*函数功能:	将输入向量各元素进行两两相减，得到输出向量
*参数：		in:传入的行向量，out传出的行向量，size传入向量的维度
*返回值：	void
*时间：		2014/11/14
**************************************************************/
void diff(const urowvec& in, irowvec& out,int size){
	for (int i = 0; i < size - 1; i++)
		out(i) = (int) (in(i + 1) - in(i));
}

/***************************************************************
*函数名:		cutST
*函数功能：	将12导联加权合成后的X or Y or Z向量进行ST-T段识别并截取
*参数:		val => X向量， sfreq => 采样频率
*			J_out => ST-T段起点(J点)，K_out => ST-T段终点（K点）
*			T_inv_flag_out => T波是否改变的标志
*返回值：	ST-T截取的段数	ST_T_count  u32类型
*时间：		2014/11/14
****************************************************************/
u32 cutST(const frowvec& val, u32 sfreq,urowvec& J_out,urowvec& K_out,bool& T_inv_flag_out)
{

	//frowvec data_in = randu<frowvec>(5);	//one of the 3 vector data
	//fvec data;
	//data.load("valX.dat", raw_ascii);

	//clock_t start = clock();

	//frowvec data_in = data.t();
	//sample frequency
	//int sfreq = 500;
	//flag indicate whether T wave is inversed
	bool Twave_inverse_flag = false;
	//smooth window
	int WW = 0;
	int p = 0;
	//threshold for recognizing R wave
	//float thresd = 0.75;
	float thresd = 0.65;		//修改thresd可以调整截取的宽度
	switch (sfreq){
	case 500:
		WW = 100;
		p = 8;
		break;
	case 1000:
		WW = 200;
		p = 16;
		break;
	default:
		Log("Error: 暂不支持的采样频率!");
		return 0;
	}
	frowvec X = val;
	//fstream fuck("D:\\1.txt", ios::out);			//测试是否X == data_in
	//(X - data_in).raw_print(fuck, "X-datain=");
	//fuck.close();
	int len = X.n_cols;
	float max_h = X.cols((int)round(len / 4), (int)round(3 * len / 4)).max();
	urowvec R_exist_area = X >= (thresd*max_h);
	//	cout << R_exist_area.cols(1,10) << endl;
	irowvec left = irowvec(R_exist_area.n_cols);
	urowvec temp = zeros<urowvec>(1);
	diff(join_horiz(temp, R_exist_area), left, R_exist_area.n_cols + 1);
	irowvec right = irowvec(R_exist_area.n_cols);
	diff(join_horiz(R_exist_area, temp), right, R_exist_area.n_cols + 1);

	//left.save("result/left.txt", raw_ascii);
	//right.save("result/right.txt", raw_ascii);

	//只能使用列向量来存储find的结果
	uvec leftind = find(left == 1);
	uvec rightind = find(right == -1);

	//leftind.save("result/leftindex.txt", raw_ascii);
	//rightind.save("result/rightindex.txt", raw_ascii);

	urowvec maxloc = zeros<urowvec>(leftind.n_rows);
	urowvec minloc = maxloc;

	for (u32 i = 0; i < leftind.n_rows; i++)
	{
		X.cols(leftind(i), rightind(i)).max(maxloc(i));
		X.cols(leftind(i), rightind(i)).min(minloc(i));
		maxloc(i) = maxloc(i) + leftind(i);
		minloc(i) = minloc(i) + leftind(i);
	}

	//求出R点的坐标集
	urowvec R_index = maxloc;
	R_index = R_index.elem(find(R_index>sfreq && R_index < (len - sfreq))).t();

	//R_index.save("result/R_index.txt", raw_ascii);
	int R_count = R_index.n_cols;
	//urowvec Ka = zeros<urowvec>(R_count - 1);
	urowvec Ka = zeros<urowvec>(R_count - 1);
	urowvec Kb = Ka, R_int = Ka;
	for (int i = 1; i < R_count; i++)
	{
		R_int(i - 1) = R_index(i) - R_index(i - 1);
		if (R_int(i - 1) < ceil(220 * sfreq / 250.0))
		{
			Ka(i-1) = R_index(i - 1) + (u32)floor(0.5*sqrt(R_int(i - 1)) + 20 * sfreq / 250.0);
			Kb(i-1) = R_index(i - 1) + (u32)floor(0.15*sqrt(R_int(i - 1)) + 30 * sfreq / 250.0);
		}
		else
		{
			Ka(i-1) = R_index(i - 1) + (u32)floor(0.5*sqrt(R_int(i - 1)) + 25 * sfreq / 250.0);
			Kb(i-1) = R_index(i - 1) + (u32)(80 * sfreq / 250.0);
		}
	}

	//Ka.save("result/Ka.txt", raw_ascii);
	//Kb.save("result/Kb.txt", raw_ascii);
	u32 RL = Ka.n_cols;
	urowvec S = zeros<urowvec>(RL);
	frowvec Ak;
	float r = 6;
	u32 akmaxind, akminind;
	for (u32 i = 0; i < RL; i++)
	{
		Ak = zeros<frowvec>(Kb(i) - Ka(i) + 1);					//原MATALB算法这里有疑问
		for (u32 k = Ka(i); k <= Kb(i); k++)
		{
			float Sk_bar = sum(X.cols(k - p, k + p)) / (2 * p + 1);
			Ak(i) = sum(X.cols(k, k + WW - 1) - Sk_bar);
		}
		Ak.max(akmaxind);
		Ak.min(akminind);
		if (abs(Ak(akminind)) > 1e-5)
		{
			float temp = abs(Ak(akmaxind)) / abs(Ak(akminind));
			if (temp> 1.0 / r && temp < r)
			{
				S(i) = max(akminind, akmaxind) + Ka(i);
			}
			else
			{
				if (abs(Ak(akminind)) > abs(Ak(akmaxind)))
					S(i) = akminind + Ka(i);
				else
					S(i) = akmaxind + Ka(i);
			}
		}
		else
		{
			S(i) = max(akminind, akmaxind) + Ka(i);

		}
		S(i)++;
	}

	//S.save("result/S.txt", raw_ascii);
	urowvec K = zeros<urowvec>(RL);
	urowvec J = K;
	for (u32 i = 0; i < RL; i++)
	{
		K(i) = S(i) + (u32)floor(180 * sfreq / 1000.0);
		//while (abs(X(K(i)))>0.005)	K(i)++;					//加上abs以准确截取双向T波
		while (X(K(i)) > 0.005)	K(i)++;
	}

	for (u32 i = 0; i < RL; i++)
	{
		int HR = (int)floor(60 * sfreq / R_int(i));
		if (HR < 100)	J(i) = R_index(i) + (u32)floor(0.6*(S(i) - R_index(i)));
		else if (HR >= 100 && HR < 110)		J(i) = R_index(i) + (u32)floor(0.55*(S(i) - R_index(i)));
		else if (110 <= HR&& HR <120)	J(i) = R_index(i) + (u32)floor(0.5*(S(i) - R_index(i)));
		else	J(i) = R_index(i) + (u32)floor(0.45*(S(i) - R_index(i)));
	}

	frowvec tempmax = zeros<frowvec>(RL);
	frowvec tempmin = tempmax;
	frowvec slope_R2J = tempmax;
	frowvec slope_J2T = tempmax;
	frowvec twave_value = tempmax;
	urowvec tempminind = zeros<urowvec>(RL);
	urowvec tempmaxind = tempminind;
	urowvec twave_index = tempminind;
	for (u32 i = 0; i < RL; i++)
	{
		tempmax(i) = X.cols(J(i), K(i)).max(tempmaxind(i));
		tempmin(i) = X.cols(J(i), K(i)).min(tempminind(i));
		if (abs(tempmax(i)) < abs(tempmin(i)))
		{
			twave_index(i) = J(i) + tempminind(i);
			twave_value(i) = tempmin(i);
		}
		else
		{
			twave_index(i) = J(i) + tempmaxind(i);
			twave_value(i) = tempmax(i);
		}
		slope_R2J(i) = (X(R_index(i)) - X(J(i))) / (R_index(i) - J(i));
		slope_J2T(i) = (X(J(i)) - twave_value(i)) / (J(i) - twave_index(i));
		if (slope_J2T(i) * slope_R2J(i) > 0)
		{
			Twave_inverse_flag = true;
			//cout << "T波倒置了!" << endl;
			break;
		}

	}

	J_out = J;
	K_out = K;
	T_inv_flag_out = Twave_inverse_flag;
	return RL;

	//clock_t end = clock();
	//cout << "所用毫秒数:" << (end - start) * 1000 / CLOCKS_PER_SEC << endl;
	//J.save("J.dat", raw_ascii);
	//K.save("K.dat", raw_ascii);
	//S.save("S.dat", raw_ascii);
	//R_index.save("R_index.dat", raw_ascii);
	//R_int.save("R_int.dat", raw_ascii);
	//Ka.save("Ka.dat", raw_ascii);
	//Kb.save("Kb.dat", raw_ascii);
	//mat R_exist_area;
	//R_exist_area.load("D:\\OurData.dat", raw_ascii);	//get the data from ascii file exported by matlab
	//cout << R_exist_area(0, 0) << endl;
	//data_in.raw_print(cout, "data_in = ");
	//cout << max_h << endl;
}


/**************************************************************
*函数名:		merge_xstt
*函数功能:	融合一维X向量截取后得到的结果
*参数:		xyz => cutST程序得到的三维向量,
*			J => J点集合, K => K点集合, x_stt_count => 截取的段数
*返回值:		x_stt => 融合之后的三维矩阵
*时间:		2015/10/4
***************************************************************/
fmat merge_xstt(const fmat& xyz, const urowvec& J, const urowvec& K,const u32& x_stt_count)
{
	int n_rows = 0;
	for (int i = 0; i < x_stt_count; i++)
	{
		n_rows += K(i) - J(i) + 1;
	}
	//cout << "merged_rows = " << n_rows << endl;
	fmat x_stt = zeros<fmat>(n_rows, 3);

	int current_block_start = 0,current_block_end = 0;
	for (int i = 0; i < x_stt_count; i++)
	{
		current_block_end = current_block_start + K(i) - J(i);
		x_stt.rows(current_block_start, current_block_end)
		 = xyz.rows(J(i), K(i));
		current_block_start = current_block_end + 1;
	}
	return x_stt;
}

/**************************************************************
*函数名:		_midf
*函数功能:	对一维列向量进行中值滤波,可配合 midfilter函数使用
*参数:		source => 待滤波的向量, result => 滤波后的向量
*			window => 滑动窗口大小
*返回值:		void
*时间:		2014/12/4
***************************************************************/
void _midf(fvec source, fvec& result,u32 window)
{
	u32 rows = source.n_rows;
	u32 n_exten, mid;
	if (window % 2 == 0)  n_exten = window / 2;
	else  n_exten = (window-1) / 2;
	fvec exten = zeros<fvec>(n_exten);
	source = join_vert(exten, source);
	source = join_vert(source, exten);
	for (u32 i = 0; i < rows; i++)
	{
		fvec temp = sort(source.rows(i, i + window - 1));
		if (window % 2 == 0)
		{
			mid = window / 2 - 1;
			result(i) = mean(temp.rows(mid, mid + 1));
		}
		else
			result(i) = temp((window - 1) / 2);
	}

}

/**************************************************************
*函数名:		midfilter1
*函数功能:	对矩阵按列向量进行中值滤波
*参数:		target => 待滤波的矩阵,window => 滤波窗口大小
*返回值:		void
*时间:		2014/12/4
***************************************************************/
void midfilter1(fmat& target,u32 window)
{
	for (u32 i = 0; i < target.n_cols; i++)
	{
		fvec result = zeros<fvec>(target.n_rows);
		_midf(target.col(i), result, window);
		target.col(i) = result;
	}
}


/**************************************************************
*函数名:		ddencmp
*函数功能:	Default values for de-noising or compression.
*参数:		dorc,worwp,signal
*返回值:		thr,sorh,keepapp
*时间:		2014/12/8
**************************************************************/
void ddencmp(string dorc, string worwp, vector<double> signal, double& thr, char& sorh, bool& keepapp)
{
	if (dorc.compare("den") == 0 && worwp.compare("wv") == 0)
	{
		//Set problem dimension.
		int dim = 1;
		//Set sorh default value.
		sorh = 's';
		// Set keepapp default value.
		keepapp = 1;

		//Set threshold default value.
		int n = signal.size();

		thr = sqrt(2 * log(n));
		vector<double> dwt_output, flag;
		vector<int> length;
		if (dim == 1)
		{
			dwt_sym(signal, 1, "db1", dwt_output, flag, length);
			dwt_output.erase(dwt_output.begin(), dwt_output.begin() + length[0]);
			for (vector<double>::iterator i = dwt_output.begin(); i != dwt_output.end(); i++)
			{
				*i = abs(*i);
			}
			sort(dwt_output.begin(), dwt_output.end());
			double normaliz = 0;
			if (dwt_output.size() % 2 != 0)
			{
				normaliz = dwt_output[(dwt_output.size() - 1) / 2];
			}
			else
			{
				normaliz = (dwt_output[(dwt_output.size() - 2) / 2] + dwt_output[dwt_output.size() / 2]) / 2;
			}

			thr = thr*normaliz / 0.6745;

		}

	}
	else
	{
		Log("invalid arguments!");
		return;
	}
}


/**************************************************************
*函数名:		wdencmp
*函数功能:	De-noising or compression using wavelets.
*参数:		o,data,len,wavename,level,flag,thr,sorh,keepapp
*返回值:		reconstructed signal
*时间:		2014/12/8
**************************************************************/
vector<double>  wdencmp(string o, vector<double>& data, vector<int> len, string wavename, int level, vector<double>& flag, double& thr, char& sorh, bool& keepapp)
{
	int dim = 1;
	int inddetstart = 0;
	int inddetend = 0;
	vector<double> cxc;
	if (o.compare("gbl") == 0)
	{
		if (keepapp)
		{
			cxc = data;
			if (dim == 1)
			{
				inddetstart = len[0];
				inddetend = data.size();
			}
			if (sorh == 's')
			{
				for (int i = inddetstart; i < data.size(); i++)
				{
					double temp = abs(data[i]) - thr;
					temp = (temp + abs(temp)) / 2;
					//cxc[i] = sign(data[i])*temp;
					if (data[i] > 0)
					{
						cxc[i] = temp;
					}
					else if (data[i] < 0) cxc[i] = -temp;
					else cxc[i] = 0;
				}
			}
		}
	}
	vector<double> xc;
	if (dim == 1)
	{
		//Wavelet reconstruction.
		idwt_sym(cxc, flag, wavename, xc, len);
	}
	return xc;
}

/**************************************************************
*函数名:		waveletfilter
*函数功能:	wavelet-denoisy & reconstruction function
*参数:		signal,level,wavename
*返回值:		reconstructed signal
*时间:		2014/12/8
**************************************************************/
vector<double> waveletfilter(vector<double> signal, int level, string wavename)
{
	vector<double> dwt_output, flag;
	vector<int> length;
	dwt_sym(signal, level, wavename, dwt_output, flag, length);
	double thr = 0;
	char sorh = 0;
	bool keepapp = 0;
	ddencmp("den", "wv", signal, thr, sorh, keepapp);
	return wdencmp("gbl", dwt_output, length, wavename, level, flag, thr, sorh, keepapp);
}

/**************************************************************
*函数名:		pretreat
*函数功能:	对原始心电信号进行预处理，生成了X、Y、Z导联，包括了中值滤波和小波变换
*参数:		filepath => 原始心电信号文件路径
*返回值:		val => 处理后的信号 sfreq => 采样频率
*时间:		2014/12/8
**************************************************************/
void pretreat(fmat& val,u32& sfreq,u32& repeat,float& TS,u32& gain)
{
	switch (val.n_cols)
	{
	case 12:
		sfreq = 500;
		val = val.rows(2000 - sfreq, 12000 + sfreq - 1);
		gain = 250;
		TS = 1.0 / sfreq;
		repeat = 5;
		break;
	default:
		Log("暂不支持的数据格式!");
	}
	fmat V3CG = zeros<fmat>(val.n_rows, 3);
	V3CG.col(0) = 0.38*val.col(0) - 0.07*val.col(1) - 0.13*val.col(6) + 0.05*val.col(7) - 0.01*val.col(8) + 0.14*val.col(9) + 0.06*val.col(10) + 0.54*val.col(11);
	V3CG.col(1) = -0.17*val.col(0) + 0.93*val.col(1) + 0.06*val.col(6) - 0.02*val.col(7) - 0.05*val.col(8) + 0.06*val.col(9) - 0.17*val.col(10) + 0.13*val.col(11);
	V3CG.col(2) = 0.11*val.col(0) + 0.23*val.col(1) - 0.43*val.col(6) - 0.06*val.col(7) - 0.14*val.col(8) - 0.20*val.col(9) - 0.11*val.col(10) + 0.31*val.col(11);
	val = join_horiz(val, V3CG);
	val = val / gain;
	fmat temp = val;
	midfilter1(temp, 200 * sfreq / 1000);
	midfilter1(temp, 600 * sfreq / 1000);
	val = val - temp;
	typedef std::vector<double> stdvec;
	for (u32 i = 0; i < 15; i++)
	{
		const stdvec& colv = conv_to< stdvec >::from(val.col(i));
		val.unsafe_col(i) = conv_to<fvec>::from(waveletfilter(colv, 7, "coif4"));
	}
	val = val.rows(sfreq, val.n_rows - sfreq - 1);
}


fmat learn(int sfreq, fmat& x_stt)
{
	clock_t start = clock();
	int n_rows = x_stt.n_rows;

	float width = 0.07;
	int num_cent = ceil(2 / width);
	int M = num_cent*num_cent*num_cent;

	fmat cent = zeros<fmat>(M, 3);

	frowvec A(1, 3);
	A.fill(1.0);

	int i, j, k, xrow = 0;
	frowvec ijk = zeros<frowvec>(1, 3);
	for (i = 1; i <= num_cent; i++)
	{
		for (j = 1; j <= num_cent; j++)
		{
			for (k = 1; k <= num_cent; k++)
			{
				xrow = (i - 1)*num_cent*num_cent + (j - 1)*num_cent + k - 1;
				ijk(0, 0) = i;
				ijk(0, 1) = j;
				ijk(0, 2) = k;
				cent.row(xrow) = A - (ijk - 1) * width;
			}
		}
	}


	int U = 3;

	frowvec x_norm = zeros<frowvec>(1, n_rows);
	for (i = 0; i < n_rows; i++)
	{
		x_norm(i) = norm(x_stt.row(i), 2);
	}
	float nor_lev = 0.9;
	float x_stt_max = x_norm.max();
	x_stt = (x_stt / x_stt_max)*nor_lev;

	float eta;
	eta = 0.8*width;

	fmat S = zeros<fmat>(M, n_rows);

	for (i = 0; i < n_rows; i++)
	{
		float m, p, q, mrow = 0;
		frowvec xyz1 = zeros<frowvec>(1, 3);
		frowvec xyz2 = zeros<frowvec>(1, 3);
		frowvec mpq = zeros<frowvec>(1, 3);
		fvec s = zeros<fvec>(M, 1);
		for (m = -U + 0.5; m <= U - 0.5; m++)
		{
			for (p = -U + 0.5; p <= U - 0.5; p++)
			{
				for (q = -U + 0.5; q <= U - 0.5; q++)
				{
					xyz1 = floor((A - x_stt.row(i)) / width) + A * 0.5;
					mpq(0, 0) = m;
					mpq(0, 1) = p;
					mpq(0, 2) = q;
					xyz2 = xyz1 + mpq;
					xyz2.elem(find(xyz2 < 1)).ones();
					xyz2.elem(find(xyz2 > num_cent)).fill(num_cent);
					mrow = (xyz2(0, 0) - 1)*num_cent*num_cent + (xyz2(0, 1) - 1)*num_cent + xyz2(0, 2) - 1;
					s(mrow, 0) = exp(-pow((norm(x_stt.row(i) - cent.row(mrow))), 2) / pow(eta, 2));
				}
			}
		}
		S.col(i) = s;

	}

	frowvec p2p = zeros<frowvec>(1, n_rows);
	frowvec p2p_sort = zeros<frowvec>(1, n_rows);
	float pp = 0;
	float p2p_30 = 0;

	for (i = 0; i < n_rows - 1; i++)
	{
		pp = norm(x_stt.row(i + 1) - x_stt.row(i));
		p2p.col(i) = pp;
	}
	p2p.col(n_rows - 1) = norm(x_stt.row(n_rows - 1) - x_stt.row(0));
	p2p_sort = arma::sort(p2p, "descend");
	p2p_30 = p2p_sort(0, 29);


	fmat W11 = zeros<fmat>(M, 1);
	fmat W21 = zeros<fmat>(M, 1);
	fmat W31 = zeros<fmat>(M, 1);
	fmat W12 = zeros<fmat>(M, 1);
	fmat W22 = zeros<fmat>(M, 1);
	fmat W32 = zeros<fmat>(M, 1);
	fmat W13 = zeros<fmat>(M, 1);
	fmat W23 = zeros<fmat>(M, 1);
	fmat W33 = zeros<fmat>(M, 1);
	mat x_hat = zeros<mat>(3, 3);

	double lambda = 10;
	double alpha = 1.99;
	double a = 0.99;
	double TS = 1.0 / sfreq;
	int repeat = 70;
	int n_rows_repeat = n_rows*repeat;

	for (i = 1; i < n_rows_repeat - 1; i++)
	{

		int ii = 0;
		int iif = 0;

		ii = i%n_rows;
		iif = ii - 1;
		switch (ii)
		{
		case 0:
			/*ii = 0;
			iif = n_rows - 1;*/
			ii = n_rows - 1;
			iif = ii - 1;
			break;
		case 1:
			iif = n_rows - 1;
		}

		if (norm(x_stt.row(ii) - x_stt.row(iif)) <= p2p_30)
		{
			double den = 0;
			den = 1 + as_scalar(lambda*S.col(iif).t()*S.col(iif));
			W12 = W11 - alpha*lambda*(x_hat(1, 0) - x_stt(ii, 0) - a*(x_hat(0, 0) - x_stt(iif, 0)))*S.col(iif) / den;
			W22 = W21 - alpha*lambda*(x_hat(1, 1) - x_stt(ii, 1) - a*(x_hat(0, 1) - x_stt(iif, 1)))*S.col(iif) / den;
			W32 = W31 - alpha*lambda*(x_hat(1, 2) - x_stt(ii, 2) - a*(x_hat(0, 2) - x_stt(iif, 2)))*S.col(iif) / den;
		}

		x_hat(2, 0) = x_stt(ii, 0) + a*(x_hat(1, 0) - x_stt(ii, 0)) + as_scalar(TS*W12.t()*S.col(ii));
		x_hat(2, 1) = x_stt(ii, 1) + a*(x_hat(1, 1) - x_stt(ii, 1)) + as_scalar(TS*W22.t()*S.col(ii));
		x_hat(2, 2) = x_stt(ii, 2) + a*(x_hat(1, 2) - x_stt(ii, 2)) + as_scalar(TS*W32.t()*S.col(ii));
		x_hat.row(0) = x_hat.row(1);
		x_hat.row(1) = x_hat.row(2);
		x_hat.row(2).fill(0);

		if (i >= (n_rows_repeat - 301))
		{
			W13 = W12 + W13;
			W23 = W22 + W23;
			W33 = W32 + W33;
		}
		W11 = W12;
		W21 = W22;
		W31 = W32;
		if (i % 10000 == 0)  {
            stringstream s;
            s << i;
            Log(s.str());
        }
	}

	fmat W = zeros<fmat>(3, M);
	fmat WS = zeros<fmat>(3, n_rows);

	W.row(0) = W13.t();
	W.row(1) = W23.t();
	W.row(2) = W33.t();

	WS = W*S / 300;

	clock_t end = clock();
    float learnTime = (end - start) / (float)CLOCKS_PER_SEC;
    stringstream s;
    s << learnTime;
	Log("learning time:" + s.str() + " seconds!");

	//WS.save("1.txt", raw_ascii);
	return WS;
}

