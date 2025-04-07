#include <iostream>
#include <vector>
#include <sstream>
#include "Rational.h"
template<int M, int N, typename Field = Rational >
class Matrix {
protected:
	std::vector<std::vector<Field>> m;
	Field determinantHelper(const std::vector<std::vector<Field>>& a, int n) const {
		if (n == 1) {
			return a[0][0];
		}
		else if (n == 2) {
			return a[0][0] * a[1][1] - a[0][1] * a[1][0];
		}
		else {
			Field det = 0;
			std::vector<std::vector<Field>> temp(n - 1, std::vector<Field>(n - 1));

			for (int p = 0; p < n; ++p) {
				int h = 0, k = 0;
				for (int i = 1; i < n; ++i) {
					for (int j = 0; j < n; ++j) {
						if (j == p) {
							continue;
						}
						temp[h][k] = a[i][j];
						k++;
						if (k == n - 1) {
							h++;
							k = 0;
						}
					}
				}
				det += a[0][p] * pow(-1, p) * determinantHelper(temp, n - 1);
			}
			return det;
		}
	}
public:
	Field det() const {
		// Ïðîâåðêà íà êâàäðàòíîñòü ìàòðèöû
		if (M != N) {
			throw std::logic_error("Determinant is only defined for square matrices.");
		}

		return determinantHelper(m, M);
	}

	Matrix() : m(M, std::vector<Field>(N, Field())) {}
	Field& operator()(size_t i, size_t j) {
		return m[i][j];
	}
	const Field& operator()(size_t i, size_t j) const{
		return m[i][j];
	}
	std::vector<Field>& operator[](size_t i){
		return m[i];
	}
	const std::vector<Field>& operator[](size_t i) const{
		return m[i];
	}
	bool operator==(const Matrix& rhs) const {
		// Ñðàâíèâàåì ðàçìåðû ìàòðèö
		if (M != rhs.m.size() || N != rhs.m[0].size()) {
			return false;
		}

		// Ïîýëåìåíòíîå ñðàâíåíèå
		for (size_t i = 0; i < M; ++i) {
			for (size_t j = 0; j < N; ++j) {
				if (m[i][j] != rhs.m[i][j]) {
					return false;
				}
			}
		}

		return true;
	}
	bool operator!=(const Matrix& rhs) const {
		return (!(*this == rhs));
	}
	Matrix& operator=(const Matrix& rhs){
		if (*this == rhs) {
			return *this;
		}
		(*this).m = rhs.m;
		return (*this);
	}
	Matrix operator-()const {
		Matrix<M, N, Field> res;
		for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				res(i, j) = -1 * (*this)(i, j);
			}
		}
		return res;
	}
	
	size_t rank() const {
		Matrix<M, N, Field> temp = *this; // Êîïèðóåì ìàòðèöó
		size_t rank = 0;

		for (size_t i = 0; i < N && rank < M; ++i) {
			// Ïîèñê ñòðîêè ñ ìàêñèìàëüíûì ýëåìåíòîì â òåêóùåì ñòîëáöå
			size_t maxRow = rank;
			for (size_t j = rank + 1; j < M; ++j) {
				if (abs(temp(j, i)) > abs(temp(maxRow, i))) {
					maxRow = j;
				}
			}

			// Åñëè ìàêñèìàëüíûé ýëåìåíò íå ðàâåí íóëþ
			if (temp(maxRow, i) != Field(0)) {
				// Ìåíÿåì ñòðîêè ìåñòàìè
				if (maxRow != rank) {
					std::swap(temp[rank], temp[maxRow]);
				}

				// Íîðìàëèçóåì òåêóùóþ ñòðîêó
				for (size_t j = i + 1; j < N; ++j) {
					temp(rank, j) /= temp(rank, i);
				}

				// Îáíóëÿåì ýëåìåíòû íèæå òåêóùåãî
				for (size_t k = 0; k < M; ++k) {
					if (k != rank && temp(k, i) != Field(0)) {
						for (size_t j = i + 1; j < N; ++j) {
							temp(k, j) -= temp(rank, j) * temp(k, i);
						}
					}
				}

				// Óâåëè÷èâàåì ðàíã
				++rank;
			}
		}

		return rank;
	}

	Matrix<M, N, Field> Mreverse() const {
		int m = M;
		if (M != N) {
			throw std::invalid_argument("Ìàòðèöà äîëæíà áûòü êâàäðàòíîé.");
		}

		// Âû÷èñëÿåì îïðåäåëèòåëü èñõîäíîé ìàòðèöû
		Field det = this->det();
		if (det == 0) {
			throw std::runtime_error("Îáðàòíîé ìàòðèöû íå ñóùåñòâóåò (îïðåäåëèòåëü ðàâåí íóëþ).");
		}

		// Ñîçäàåì ðåçóëüòèðóþùóþ ìàòðèöó
		Matrix<M, N, Field> rez;

		// Çàïîëíÿåì ðåçóëüòèðóþùóþ ìàòðèöó àëãåáðàè÷åñêèìè äîïîëíåíèÿìè
		for (size_t i = 0; i < m; i++) {
			for (size_t j = 0; j < m; j++) {
				// Âû÷èñëÿåì ìèíîð äëÿ ýëåìåíòà (i, j)
				std::vector<std::vector<Field>> minorMat = GetMatr(this->m, i, j);
				Field minor = determinantHelper(minorMat, m - 1);

				// Âû÷èñëÿåì àëãåáðàè÷åñêîå äîïîëíåíèå
				if ((i + j) % 2 == 1) {
					minor =  minor * (-1);  // Ìåíÿåì çíàê, åñëè ñóììà èíäåêñîâ íå÷¸òíàÿ
				}

				// Çàïèñûâàåì ýëåìåíò îáðàòíîé ìàòðèöû (ñðàçó òðàíñïîíèðóåì)
				rez(j, i) = minor / det;
			}
		}

		return rez;
	}

	Matrix& operator+=(const Matrix rhs) {
		for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				(*this)(i, j) += rhs(i, j);
			}
		}
		return *this;
	}
	Matrix& operator-=(const Matrix& rhs) {
		for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				(*this)(i, j) -= rhs(i, j);
			}
		}
		return *this;
	}
	friend std::ostream& operator<<(std::ostream& os, const Matrix<M, N, Field>& a) {
		for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				os << a(i, j) << " ";
			}
			os << std::endl;
		}
		return os;
	}
	friend std::istream& operator>>(std::istream& in, Matrix<M, N, Field>& matr) {
		
		for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				in >> matr[i][j];
			}
		}
		return in;
	}
	friend Matrix<M, N, Field> operator-(Matrix<M, N, Field>& a, Matrix<M, N, Field>& b) {
		Matrix<M, N, Field> c;
		for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				c(i, j) = a(i, j) - b(i, j);
			}
		}
		return c;
	}
	friend Matrix<M, N, Field> operator+(Matrix<M, N, Field>& a, Matrix<M, N, Field>& b) {
		Matrix<M, N, Field> c;
		for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				c(i, j) = a(i, j) + b(i, j);
			}
		}
		return c;
	}
	
	Matrix<N, M, Field> transposed() const{
		Matrix<N, M, Field> a;
		for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				a(j, i) = (*this)(i, j);
			}
		}
		
		return a;
	}
	Matrix<M, N, Field>& operator*(Field a) {
		for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				(*this)(i, j) *= a;
			}
		}
		return *this;
	}
	friend Matrix<M, N, Field>& operator*(Field a, Matrix<M, N, Field>& rhs) {
		for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				rhs(i, j) *= a;
			}
		}
		return rhs;
	}
	std::vector<Field>& GetRows(size_t ia) {
		return (*this)[ia];
	}
	std::vector<Field> GetCols(size_t ja) {
		std::vector<Field> col;
		for (size_t i = 0; i < M; i++)
		{
			col.push_back((*this)[i][ja]);
		}
		return col;
	}
	Field trace(const Matrix<N, N, Field>& mat) {

		Field sum = 0;
		for (size_t i = 0; i < N; i++)
		{
			sum += m[i][i];
		}
		return sum;
	}

	

	std::vector<std::vector<Field>> GetMatr(const std::vector<std::vector<Field>>& mat, size_t excludeRow, size_t excludeCol) const {
		size_t siz = mat.size();
		std::vector<std::vector<Field>> minorMat(siz - 1, std::vector<Field>(siz - 1));

		size_t minorRow = 0, minorCol = 0;
		for (size_t i = 0; i < siz; i++) {
			if (i == excludeRow) continue;
			for (size_t j = 0; j < siz; j++) {
				if (j == excludeCol) continue;
				minorMat[minorRow][minorCol] = mat[i][j];
				minorCol++;
			}
			minorRow++;
			minorCol = 0;
		}

		return minorMat;
	}
};

template<int M, int N, int K, typename Field>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field>& a, const Matrix<N, K, Field>& b) {
	Matrix<M, K, Field> c;
	for (size_t i = 0; i < M; i++)
	{
		for (size_t k = 0; k < K; k++) {
			c(i, k) = Field();
			for (size_t j = 0; j < N; j++)
			{
				c[i][k] += a[i][j] * b[j][k];
			}
		}


	}
	return c;
}
	





template<size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;



int main() {
	SquareMatrix<3> mat;
	mat(0, 0) = 1; mat(0, 1) = 2; mat(0, 2) = 3;
	mat(1, 0) = 4; mat(1, 1) = 5; mat(1, 2) = 6;
	mat(2, 0) = 7; mat(2, 1) = 8; mat(2, 2) = 9;

	SquareMatrix<3> mat1;
	mat1(0, 0) = 2; mat1(0, 1) = 8; mat1(0, 2) = 4;
	mat1(1, 0) = 5; mat1(1, 1) = 1; mat1(1, 2) = 7;
	mat1(2, 0) = 1; mat1(2, 1) = 2; mat1(2, 2) = 9;

	std::cout << "Matrix:\n" << mat << std::endl;
	std::cout << "Matrix1:\n" << mat1 << std::endl;
	std::cout << "Determinant: " << mat.det() << std::endl;
	std::cout << "Rank: " << mat.rank() << std::endl;
	std::cout << mat * 2 << std::endl;

	if (mat1 != mat) {
		std::cout << "YES" << std::endl;
	}

	std::cout << "Inverse matrix:\n" << mat1.Mreverse() << std::endl;

	return 0;
}
