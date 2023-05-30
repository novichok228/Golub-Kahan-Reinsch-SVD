#include <iostream>
#include <fstream>
#include <typeinfo>
#include <string>
#include <ctime>
#include <math.h>
#include <iomanip>
#include <vector>
#include<time.h>
#include <algorithm>
#include <random>
#include <iostream>
#include <stack>
#include <stdlib.h>
constexpr double eps = 0.000001;// точность вычисления сингулярных значений и сравнения элементов матриц
using namespace std;
class Exception : public std::exception
{
protected:
	char* str;
public:
	Exception(const char* s)
	{
		str = new char[strlen(s) + 1];
		strcpy_s(str, strlen(s) + 1, s);
	}
	Exception(const Exception& e)
	{
		str = new char[strlen(e.str) + 1];
		strcpy_s(str, strlen(e.str) + 1, e.str);
	}
	virtual ~Exception()
	{
		if (str != NULL) { delete[] str; }
	}
	virtual void print() const
	{
		cout << "\nException: " << str;
	}
};

class IndexOutOfBoundsException : public Exception {
public:
	IndexOutOfBoundsException(const char* s) : Exception(s) {}
	IndexOutOfBoundsException(const IndexOutOfBoundsException& e) : Exception(e) {}
	void print() const {
		cout << "\nIndexOutOfBoundsException: " << str;
	}
	virtual ~IndexOutOfBoundsException() {}
};

class WrongDimensionException : public Exception {
public:
	WrongDimensionException(const char* s) : Exception(s) {}
	WrongDimensionException(const WrongDimensionException& e) : Exception(e) {}
	void print() const {
		cout << "\n WrongDimensionException: " << str;
	}
	virtual ~WrongDimensionException() {}
};

class WrongSizeException : public WrongDimensionException {
public:
	WrongSizeException(const char* s) : WrongDimensionException(s) {}
	WrongSizeException(const WrongSizeException& e) : WrongDimensionException(e) {}
	void print() const {
		cout << "\n WrongSizeException: " << str;
	}
	virtual ~WrongSizeException() {}
};

class NotEnoughSpaceException : public Exception {
public:
	NotEnoughSpaceException(const char* s) : Exception(s) {}
	NotEnoughSpaceException(const NotEnoughSpaceException& e) : Exception(e) {}
	void print() const {
		cout << "\n NotEnoughSpaceException" << str;
	}
};

template<class Type>
class BaseMatrix
{
protected:
	Type** ptr;
	int height;
	int width;
public:
	BaseMatrix<Type>(int Height = 2, int Width = 2)
	{
		if (Height <= 0 || Width <= 0) { throw WrongSizeException("BaseMatrixConstructor: Non-positive size of matrix"); }
		height = Height;
		width = Width;
		if (ptr = new Type * [height]) {
			for (int i = 0; i < height; ++i) {
				if (!(ptr[i] = new Type[width])) { throw NotEnoughSpaceException("BaseMatrixConstructor: NotEnoughSpace 1"); }
			}
		}
		else { throw NotEnoughSpaceException("BaseMatrixConstructor: NotEnoughSpace 2"); }
		
		for (int i = 0; i < height; ++i)//создаем нулевую матрицу
			for (int j = 0; j < width; ++j)
				ptr[i][j] = 0;
	}
	BaseMatrix<Type>(const BaseMatrix<Type>& M)
	{
		this->height = M.height;
		this->width = M.width;
		if (ptr = new Type * [height]) {
			for (int i = 0; i < height; ++i)
			{
				if (ptr[i] = new Type[width]) {
					for (int j = 0; j < width; ++j) { ptr[i][j] = M.ptr[i][j]; }
				}
				else { throw NotEnoughSpaceException("BaseMatrixCopyConstructor: NotEnoughSpace 1"); }
			}
		}
		else { throw NotEnoughSpaceException("BaseMatrixCopyConstructor: NotEnoughSpace 2"); }
	}
	void I() {//пологаем главную диагональ еденицами(предполагаем, что эта функция создает единичную матрицу, так как работаем изначально с нулевыми)
		int len = min(height, width);
		for (int i = 0; i < len; ++i)
			ptr[i][i] = 1;
	}
	BaseMatrix<Type> operator+(const BaseMatrix<Type>& M) {
		BaseMatrix<Type> A(this->height, this->width);
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				A.ptr[i][j] = ptr[i][j] + M.ptr[i][j];
			}
		}
		return A;
	}
	bool operator==(const BaseMatrix<Type>& M) const{
		if (height != M.height || width != M.width) throw WrongSizeException("operator==: wrong sizes");
		Type fraben_norm1 = 0;
        Type fraben_norm2 = 0;
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                fraben_norm1 += (ptr[i][j] * ptr[i][j]);
                fraben_norm2 += (M.ptr[i][j] * M.ptr[i][j]);
            }
        }
        fraben_norm1 = sqrt(fraben_norm1);
        fraben_norm2 = sqrt(fraben_norm2);
        if (abs(fraben_norm1 - fraben_norm2) >= eps) return false;
		return true;
	}
	BaseMatrix<Type> operator-(const BaseMatrix<Type>& M) const{
		BaseMatrix<Type> A(this->height, this->width);
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				A.ptr[i][j] = ptr[i][j] - M.ptr[i][j];
			}
		}
		return A;
	}
	Type mult(const BaseMatrix<Type>& M) const {//скалярное произведение матриц -- возвращает число
		if (M.height > 1 && M.width > 1) { throw WrongDimensionException("mult1: wrong dimension"); }
		if (this->height != M.height || this->width != M.width) { throw WrongDimensionException("Mult2"); }
		Type res = 0;
		for (int i = 0; i < this->height; ++i) {
			for (int j = 0; j < this->width; ++j) {
				res += this->ptr[i][j] * M.ptr[i][j];
			}
		}
		return res;
	}
	BaseMatrix<Type>& operator=(const BaseMatrix<Type>& M) {
		if (this == &M) return *this;
		if (ptr != NULL) {
			for (int i = 0; i < height; ++i) {
				if (ptr[i] != NULL) {
					delete[] ptr[i];
				}
			}
			delete[] ptr;
		}
		this->height = M.height;
		this->width = M.width;
		if (ptr = new Type * [height]) {
			for (int i = 0; i < height; ++i) {
				if (ptr[i] = new Type[width]) {
					for (int j = 0; j < width; ++j) {
						ptr[i][j] = M.ptr[i][j];
					}
				}
				else { throw NotEnoughSpaceException("BaseMatrixOperator=: NotEnoughSpace 1"); }
			}
		}
		else { throw NotEnoughSpaceException("BaseMatrixOperator=: NotEnoughSpace 2"); }
		return *this;
	}
	Type* operator[](int index)
	{
		if (index > height) { throw WrongSizeException("operator[1]: wrong size"); }
		if (index < 0) { throw WrongSizeException("operator[2]: wrongsize"); }
		return ptr[index];
	}

	Type& operator()(int index1, int index2) {
		if (index1 > height || index2 > width)
		{
			throw WrongSizeException("operator(1): wrongsize");
		}
		if (index1 < 0 || index2 < 0)
		{
			throw WrongSizeException("operator(2): wrongsize");
		}

		return ptr[index1][index2];
	}
	BaseMatrix<Type> operator*(const BaseMatrix<Type>& M) const {
		if (this->width != M.height) //{ throw WrongSizeException("operator*: wrong sizes"); }
		{//special for householder decomposition
			if (M.height < this->width) {
				BaseMatrix<Type> M2(this->width, M.width);
				for (int i = 0; i < M.height; ++i) {
					for (int j = 0; j < M.width; ++j) {
						M2.ptr[i][j] = M.ptr[i][j];
					}
				}
				BaseMatrix<Type> C(this->height, M2.width);

				for (int i = 0; i < this->height; ++i) {
					for (int j = 0; j < M2.width; ++j) {
						for (int k = 0; k < M2.height; ++k)
							C[i][j] += this->ptr[i][k] * M2.ptr[k][j];
					}
				}
				return C;
			}
		}

		BaseMatrix<Type> C(this->height, M.width);

		for (int i = 0; i < this->height; ++i) {
			for (int j = 0; j < M.width; ++j) {
				for (int k = 0; k < M.height; ++k)
					C[i][j] += this->ptr[i][k] * M.ptr[k][j];
			}
		}
		return C;
	}
	BaseMatrix<Type>& operator*(const Type c) {
		for (int i = 0; i < this->height; ++i) {
			for (int j = 0; j < this->width; ++j) {
				this->ptr[i][j] *= c;
			}
		}
		return *this;
	}
	BaseMatrix<Type>& operator/(const Type c) {
		for (int i = 0; i < this->height; ++i) {
			for (int j = 0; j < this->width; ++j) {
				this->ptr[i][j] /= c;
			}
		}
		return *this;
	}
	int get_height() const {
		return height;
	}
	int get_width() const {
		return width;
	}
	void zero() {//положить все элементы матрицы нулями
		for (int i = 0; i < this->height; ++i)
			for (int j = 0; j < this->width; ++j) {
				this->ptr[i][j] = 0;
			}
	}

	BaseMatrix<Type>& round() {//округлить околонулевые элементы
		for (int i = 0; i < this->height; i++) {
			for (int j = 0; j < this->width; j++) {
				ptr[i][j] = (abs(ptr[i][j]) < eps) ? 0 : ptr[i][j];
			}
		}
		return *this;
	}
	BaseMatrix<Type> column(const int j) const {//получить j-ю колонку
		BaseMatrix<Type> res(height, 1);
		for (int i = 0; i < height; ++i) {
			res[i][0] = ptr[i][j];
		}
		return res;
	}
	BaseMatrix<Type> column2(int start, int finish, const int j) const {//получить j-ю колонку от start элемента до finisf
		BaseMatrix<Type> res(height, 1);
		for (int i = start; i < height-finish; ++i) {
			res[i][0] = ptr[i][j];
		}
		return res;
	}
	BaseMatrix<Type> Tm()const {//транспонировать матрицу
		BaseMatrix<Type> res(this->width, this->height);

		for (int i = 0; i < this->height; i++) {
			for (int j = 0; j < this->width; j++) {
				res.ptr[j][i] = this->ptr[i][j];
			}
		}
		return res;
	}
	bool is_super_diag() const {//Проверка на то, что матрица является верхне-бидиагональной
		for (int i = 0; i < height - 1; ++i) {
			for (int j = 0; j < width; ++j) {
				if ((i == j) || (j == i + 1)) continue;
				else if (ptr[i][j] != 0) {
					cout << i << " " << j;
					return false;
				}
			}
		}
		return true;
	}
	bool get_diag_is_diag(int k) const//работает с конца || только квадратные бидиагональные матрицы
	{
		
		int i = 1;
		while ((i < k) && (ptr[height - i - 1][height - i] != 0)) {
			i++;
			return false;
		}
		return true;
	}
	bool get_diag_is_superdiag(int q, int k) const//работает с q-го элемента || только квадратные бидиагональные матрицы
	{

		int i = 1;
		while ((i < k) && (ptr[height - i - 1-q][height - i-q] == 0)) {
			i++;
			return false;
		}
		return true;
	}
	BaseMatrix<Type> clear_firsts(int k) const//удалить первые к колонок и рядов
	{
		BaseMatrix<Type> res(height - k, width - k);
		for (int i = k; i < height; ++i) {
			for (int j = k; j < width; ++j) {
				res[i - k][j - k] = ptr[i][j];
			}
		}
		return res;
	}
	BaseMatrix<Type> clear_last_row(int k) const//удалить последние к рядов
	{
		BaseMatrix<Type> res(height - k, width);
		for (int i = 0; i < height - k; ++i) {
			for (int j = 0; j < width; ++j) {
				res[i][j] = ptr[i][j];
			}
		}
		return res;
	}
	BaseMatrix<Type> clear_lasts(int k) const {//удалить последние к колонок и рядов
		BaseMatrix<Type> res(height - k, width - k);
		for (int i = 0; i < height - k; ++i) {
			for (int j = 0; j < width - k; ++j) {
				res[i][j] = ptr[i][j];
			}
		}
		return res;
	}
	BaseMatrix<Type> expand_last_row(int k) const {//расширить матрицу на к рядов
		BaseMatrix<Type> res(height + k, width);
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				res.ptr[i][j] = ptr[i][j];
			}
		}
		for (int i = height; i < height + k; ++i) {
			for (int j = 0; j < width; ++j) {
				res.ptr[i][j] = 0;
			}
		}
		return res;
	}
	friend ostream& operator<<(ostream& ustream, const BaseMatrix<Type> obj)
	{
		if (typeid(ustream).name() == typeid(ofstream).name())
		{

			ustream << obj.height << " " << obj.width << "\n";
			for (int i = 0; i < obj.height; i++)
			{
				for (int j = 0; j < obj.width; j++)
					ustream << obj.ptr[i][j] << "                ";
				ustream << "\n";
			}

			return ustream;
		}
		for (int i = 0; i < obj.height; i++)
		{
			for (int j = 0; j < obj.width; j++)
				ustream << obj.ptr[i][j] << "       ";
			ustream << "\n";
		}
		return ustream;
	}
	virtual ~BaseMatrix<Type>()
	{
		if (ptr != NULL)
		{
			for (int i = 0; i < height; ++i) {
				if (ptr[i] != NULL) { delete[] ptr[i]; }
			}
			delete[] ptr;
		}
	}
	void print() const
	{
		cout << height << ' ' << width << '\n';
	}
};


template <class t>
void householderBidiagonalization(BaseMatrix<t>& Q, BaseMatrix<t>& R, BaseMatrix<t>& S) {
	t mag, alpha;
	int m = R.get_height(), n = R.get_width();
	BaseMatrix<t> u(m, 1), v(m, 1),
		u_(n, 1), v_(n, 1);

	BaseMatrix<t> P(m, m), I(m, m),
		P_(n, n), I_(n, n);
	P.I(); I.I(); P_.I(); I_.I();
	Q = BaseMatrix<t>(m, m);
	Q.I();
	S = BaseMatrix<t>(n, n);
	S.I();

	for (int i = 0; i < n; i++) {
		u.zero(); v.zero();

		mag = 0.0;
		for (int j = i; j < m; j++) {
			u[j][0] = R[j][i];
			mag += u[j][0] * u[j][0];
		}
		mag = sqrt(mag);

		alpha = u[i][0] < 0 ? mag : -mag;

		mag = 0.0;
		for (int j = i; j < m; j++) {
			v[j][0] = j == i ? u[j][0] + alpha : u[j][0];
			mag += v[j][0] * v[j][0];
		}
		mag = sqrt(mag);

		if (mag > eps) {
			for (int j = i; j < m; j++) v[j][0] /= mag;
			R = R - v * (v.Tm() * R) * 2.0;// (mx1) * ((1xm) *(m*n)) \\ (mx1)*(1xn)  m*m*n+n
			Q = Q - (Q * v) * v.Tm() * 2.0;
		}
		u_.zero(); v_.zero();

		mag = 0.0;
		for (int j = i + 1; j < n; j++) {
			u_[j][0] = R[i][j];
			mag += u_[j][0] * u_[j][0];
		}

		mag = sqrt(mag);

		if ((i != n - 1)) alpha = u_[i + 1][0] < 0 ? mag : -mag;
		else alpha = mag;
		mag = 0.0;
		for (int j = i + 1; j < n; j++) {
			v_[j][0] = j == i + 1 ? u_[j][0] + alpha : u_[j][0];
			mag += v_[j][0] * v_[j][0];
		}
		mag = sqrt(mag);

		if (mag > eps) {
			for (int j = i + 1; j < n; j++) v_[j][0] /= mag;
			R = R - (R * v_) * v_.Tm() * 2.0;
			S = S - (S*v_ )* v_.Tm() * 2.0;// (nxn * nx1) * 1xn
		}
	}
}


template <class t>
void delete_zero(BaseMatrix<t>& B, int row) {//если у бидиагональной матрицы на главной диагонали встречается 0 -- эта функция удаляет 0 на всем ряду
	int n = B.get_height();
	t alpha, beta, c, s, r1, r2;
	for (int i = 0; i < n - row - 1; ++i) {
		alpha = B[row][row+1+i];
		beta = B[row+1+i][row+1+i];
		c = beta / sqrt(alpha * alpha + beta * beta);
		s = alpha / sqrt(alpha * alpha + beta * beta);
		for (int j = 0; j < n; ++j) {
			r1 = B[row][j];
			r2 = B[row+ 1+i][j];
			B[row][j] = c * r1 - s * r2;
			B[row+i + 1][j] = s * r1 + c * r2;
		}
	}
}

template <class t>
void Golub_Reinsch(BaseMatrix<t>& Q1, BaseMatrix<t>& B, BaseMatrix<t>& S) {
	int m = B.get_height();
	int n = B.get_width();
	B = B.clear_last_row(m - n);
	m = n;
	while (true) {
		for (int i = 0; i < n - 1; ++i) {
			if ((abs(B[i][i + 1])) <= eps * (abs(B[i][i]) + abs(B[i + 1][i + 1]))) {
				B[i][i + 1] = 0;
			}
		}
		
		int n_p_q = 0;
		
		int q = 2;
		while ((q <= n) && (B.get_diag_is_diag(q))) {
			q++;
		}
		q--;
		if (q == n) {
			break;
		}
		while ((q != 0) && (B[n - q - 1][n - q] != 0)) {
			q--;
		}
		
		int i = 2;
		while ((i < n- q +1)&&(B.get_diag_is_superdiag(q,i))) {
			n_p_q = i;
			i++;
		}
		int p = n - n_p_q - q;
		bool flag = false;
		for (int i = p; i < n - q - 1; ++i) {
			if (B[i][i] == 0) {
				delete_zero(B, i);
				flag = true;
			}
		}
		if (flag == false && n_p_q != 0) {
			Golub_kahan(Q1, B, S, p, q, n);
		}

	}
	m = Q1.get_width();
	B = B.expand_last_row(m - n);
}

template <class t>
void Golub_kahan(BaseMatrix<t>& Q1, BaseMatrix<t>& B, BaseMatrix<t>& S, int p, int q, int n) {
	int m = Q1.get_width();
	BaseMatrix<t> col1 = B.column2(p,q,B.get_height() - 2);
	BaseMatrix<t> col2 = B.column2(p,q,B.get_width() - 1);
	t c00 = col1.mult(col1);
	t c01 = col1.mult(col2);
	t c10 = col2.mult(col1);
	t c11 = col2.mult(col2);
	t u = c11;
	t u1 = (c11 + c00 - sqrt((c11 + c00) * (c11 + c00) - 4 * (c00 * c11 - c01 * c10))) / 2;
	t u2 = (c11 + c00 + sqrt((c11 + c00) * (c11 + c00) - 4 * (c00 * c11 - c01 * c10))) / 2;
	u = (abs(u1 - u) < abs(u2 - u)) ? u1 : u2;
	t alpha = B[p][p] * B[p][p] - u;
	t beta = B[p][p] * B[p][p + 1];
	t r1, r2, c, s;
	for (int k = p; k < n - q - 1; ++k) {
		c = alpha / sqrt(alpha * alpha + beta * beta);
		s = -beta / sqrt(alpha * alpha + beta * beta);
		for (int j = 0; j < n; ++j) {
			r1 = B[j][k];
			r2 = B[j][k + 1];
			B[j][k] = c * r1 - s * r2;
			B[j][k + 1] = s * r1 + c * r2;
			r1 = S[j][k];
			r2 = S[j][k + 1];
			S[j][k] = c * r1 - s * r2;
			S[j][k + 1] = s * r1 + c * r2;
		}
		alpha = B[k][k];
		beta = B[k + 1][k];
		c = alpha / sqrt(alpha * alpha + beta * beta);
		s = -beta / sqrt(alpha * alpha + beta * beta);
		for (int j = 0; j < n; ++j) {
			r1 = B[k][j];
			r2 = B[k + 1][j];
			B[k][j] = c * r1 - s * r2;
			B[k + 1][j] = s * r1 + c * r2;
		}
		for (int j = 0; j < m; ++j) {
			r1 = Q1[j][k];
			r2 = Q1[j][k + 1];
			Q1[j][k] = c * r1 - s * r2;
			Q1[j][k + 1] = s * r1 + c * r2;
		}
		if (k < n - q - 1) {
			alpha = B[k][k + 1];
			beta = B[k][k + 2];
		}
	}

}
template <class t>
void SVD(BaseMatrix<t>& a) {
	BaseMatrix<double> Q1;
	BaseMatrix<double> R1;
	BaseMatrix<double> S;
	R1 = a;

	srand(time(NULL));
	float timer3 = std::clock();
	cout << "start\n";

	householderBidiagonalization(Q1, R1, S);
	std::cout << std::endl << "bidiagonalization time: " << (float)(clock() - timer3) / CLOCKS_PER_SEC << std::endl;
	Golub_Reinsch(Q1, R1, S);
	std::cout << std::endl << "svd completed(bidiagonalize + diagonalize) time: " << (float)(clock() - timer3) / CLOCKS_PER_SEC << std::endl;


	cout << "\n\n" << ((Q1 * R1 * S.Tm()).round() == a) << "\n\n";//проверка, что Q*R*S' == A
	BaseMatrix<double> Id(Q1.get_height(), Q1.get_width());
	BaseMatrix<double> Id2(S.get_height(), S.get_width());
	Id.I();
	Id2.I();


	BaseMatrix<double> Q = (Q1 * Q1.Tm()).round();//проверка на унитарность Q
	cout << (Q == Id) << "\n\n";
	Q = (Q1.Tm() * Q1).round();
	cout << (Q == Id) << "\n\n";
	BaseMatrix<double> S2 = (S * S.Tm()).round();//проверка на унитарность S
	cout << (S2 == Id2) << "\n\n";
	S2 = (S.Tm() * S).round();
	cout << (S2 == Id2) << "\n\n";
	BaseMatrix<double> R = R1.round();//проверка на двухдиагональность R
	cout << R1.is_super_diag();


	//cout << (Q1.Tm() * Q1).round() << "\n\n";
	//cout << (S.Tm() * S).round() << "\n\n";
	//cout << (S * S.Tm()).round() << "\n\n";
	//cout << (R1).round() << "\n\n";
}

int main()
{
	int m = 3000, n = 2000; //set matrix size
	const int RANDVALUE = 1000;
	BaseMatrix<double> a(m, n);
	
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			a[i][j] = (rand() % RANDVALUE - RANDVALUE / 2);
		}
	}
	SVD(a);
	
}
