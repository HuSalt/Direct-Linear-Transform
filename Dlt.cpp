#include "Dlt.h"

//constructor
Dlt::Dlt()
{
	//�I���������Y�]�w
	x = new double[4];
	y = new double[4];
	xt = new double[4];
	yt = new double[4];
}

Dlt::Dlt(double * in_x, double * in_y, double * in_xt, double * in_yt)
{
	//�I���������Y�]�w
	x = new double[4];
	y = new double[4];
	xt = new double[4];
	yt = new double[4];

	for (int i = 0; i < 4; i++)
	{
		x[i] = in_x[i];
		y[i] = in_y[i];
		xt[i] = in_xt[i];
		yt[i] = in_yt[i];
	}
}

//destructor
Dlt::~Dlt()
{
	delete[] x;
	delete[] y;
	delete[] xt;
	delete[] yt;
}


//�^�ǥ��W�ƹL��Dlt
void Dlt::normalize(Mat & o_T, Mat & t_T, Dlt& normal)
{
	//���W�Ưx�}�ѼƷǳ�
	double x_av = 0, y_av = 0, o_norm = 0, o_norm_y = 0;
	double xt_av = 0, yt_av = 0, t_norm = 0, t_norm_y = 0;

	//�`�M>>����
	for (int i = 0; i < 4; i++)
	{
		x_av += x[i];
		y_av += y[i];

		xt_av += xt[i];
		yt_av += yt[i];
	}

	//�U�I������>>�����I
	x_av /= 4;
	y_av /= 4;

	xt_av /= 4;
	yt_av /= 4;

	//�Z���M>>����
	for (int i = 0; i < 4; i++)
	{
		double tmp;
		tmp = abs(x[i] - x_av);
		o_norm += tmp;

		tmp = abs(y[i] - y_av);
		o_norm_y += tmp;

		tmp = abs(xt[i] - xt_av);
		t_norm += tmp;

		tmp = abs(yt[i] - yt_av);
		t_norm_y += tmp;
	}

	//���O��x,y���W��(mapping ���@���ץ����Z��(x,y)��1 >> �G���U�����Z���ڸ�2)
	o_norm /= 4;
	t_norm /= 4;
	o_norm_y /= 4;
	t_norm_y /= 4;

	//���W�Ưx�}
	o_T = (Mat_<double>(3, 3) <<
		o_norm, 0, -o_norm * x_av,
		0, o_norm_y, -o_norm_y * y_av,
		0, 0, 1);

	t_T = (Mat_<double>(3, 3) <<
		t_norm, 0, -t_norm * xt_av,
		0, t_norm_y, -t_norm_y * yt_av,
		0, 0, 1);

	//�t�~�x�s���W�ƫ᪺DLT���k
	double* n_x = new double[4];
	double* n_y = new double[4];
	double* n_xt = new double[4];
	double* n_yt = new double[4];

	//���W�Ʃһ��I
	for (int i = 0; i < 4; i++)
	{
		//�I�s�Jmat
		Mat pt = (Mat_<double>(3, 1) <<
			x[i], y[i], 1);
		//���W��A4�ȤW���I
		pt = o_T * pt;


		n_x[i] = pt.at<double>(0, 0);
		n_y[i] = pt.at<double>(1, 0);


		//�I�s�Jmat
		pt = (Mat_<double>(3, 1) <<
			xt[i], yt[i], 1);
		//���W�ƬM�g���I
		pt = t_T * pt;

		n_xt[i] = pt.at<double>(0, 0);
		n_yt[i] = pt.at<double>(1, 0);
	}

	//Dlt temp(n_x, n_y, n_xt, n_yt);
	normal.Dlt::Dlt(n_x, n_y, n_xt, n_yt);

	//���ϯx�}�A�Ω󥿳W���٭�L�{
	invert(t_T, t_T);
}

//��l��DLT�һݪ�A�x�}
Mat Dlt::mat_Init()
{
	//DLT���x�}A
	Mat A = (Mat_<double>(12, 9) <<
		0, 0, 0, -getValue('x',false,0), -getValue('y', false, 0), -1, getValue('y', true, 0) * getValue('x', false, 0), getValue('y', true, 0) * getValue('y', false, 0), getValue('y', true, 0),
		getValue('x', false, 0), getValue('y', false, 0), 1, 0, 0, 0, -getValue('x', true, 0) * getValue('x', false, 0), -getValue('x', true, 0) * getValue('y', false, 0), -getValue('x', true, 0),
		-getValue('y', true, 0) * getValue('x', false, 0), -getValue('y', true, 0) * getValue('y', false, 0), -getValue('y', true, 0), getValue('x', true, 0) * getValue('x', false, 0), getValue('x', true, 0) * getValue('y', false, 0), getValue('x', true, 0), 0, 0, 0,

		0, 0, 0, -getValue('x', false, 1), -getValue('y', false, 1), -1, getValue('y', true, 1) * getValue('x', false, 1), getValue('y', true, 1) * getValue('y', false, 1), getValue('y', true, 1),
		getValue('x', false, 1), getValue('y', false, 1), 1, 0, 0, 0, -getValue('x', true, 1) * getValue('x', false, 1), -getValue('x', true, 1) * getValue('y', false, 1), -getValue('x', true, 1),
		-getValue('y', true, 1) * getValue('x', false, 1), -getValue('y', true, 1) * getValue('y', false, 1), -getValue('y', true, 1), getValue('x', true, 1) * getValue('x', false, 1), getValue('x', true, 1) * getValue('y', false, 1), getValue('x', true, 1), 0, 0, 0,

		0, 0, 0, -getValue('x', false, 2), -getValue('y', false, 2), -1, getValue('y', true, 2) * getValue('x', false, 2), getValue('y', true, 2) * getValue('y', false, 2), getValue('y', true, 2),
		getValue('x', false, 2), getValue('y', false, 2), 1, 0, 0, 0, -getValue('x', true, 2) * getValue('x', false, 2), -getValue('x', true, 2) * getValue('y', false, 2), -getValue('x', true, 2),
		-getValue('y', true, 2) * getValue('x', false, 2), -getValue('y', true, 2) * getValue('y', false, 2), -getValue('y', true, 2), getValue('x', true, 2) * getValue('x', false, 2), getValue('x', true, 2) * getValue('y', false, 2), getValue('x', true, 2), 0, 0, 0,


		0, 0, 0, -getValue('x', false, 3), -getValue('y', false, 3), -1, getValue('y', true, 3) * getValue('x', false, 3), getValue('y', true, 3) * getValue('y', false, 3), getValue('y', true, 3),
		getValue('x', false, 3), getValue('y', false, 3), 1, 0, 0, 0, -getValue('x', true, 3) * getValue('x', false, 3), -getValue('x', true, 3) * getValue('y', false, 3), -getValue('x', true, 3),
		-getValue('y', true, 3) * getValue('x', false, 3), -getValue('y', true, 3) * getValue('y', false, 3), -getValue('y', true, 3), getValue('x', true, 3) * getValue('x', false, 3), getValue('x', true, 3) * getValue('y', false, 3), getValue('x', true, 3), 0, 0, 0
		);
	
	return A;
}

//�p��H�x�}
Mat Dlt::getMatH()
{
	Mat A, H;
	A = mat_Init();

	//�DH�x�}(��A��SVD�o��M�g�x�}H)
	Mat w, u, v;
	SVD::compute(A, w, u, v);

	H = (Mat_<double>(3, 3) <<
		v.at<double>(8, 0), v.at<double>(8, 1), v.at<double>(8, 2),
		v.at<double>(8, 3), v.at<double>(8, 4), v.at<double>(8, 5),
		v.at<double>(8, 6), v.at<double>(8, 7), v.at<double>(8, 8)
		);

	A.release();
	w.release();
	u.release();
	v.release();

	
	return H;
}

//���Ȩ��
double Dlt::getValue(char c, bool trans, int index)
{
	if (c == 'x')
	{
		if (trans)
		{
			return xt[index];
		}
		return x[index];
	}
	else if (c == 'y')
	{
		if (trans)
		{
			return yt[index];
		}
		return y[index];
	}
	else
	{
		return 0.0;
	}
}