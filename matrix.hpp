#ifndef MATRIX_HPP
#define MATRIX_HPP

namespace caams
{
	enum matrix_flags{
		uninitialized = 0,
		init = 1,
		identity = 2,
		init_identity = 3,
		zeros = 1
	};
	struct matrix
	{
	public:
		long rows;
		long cols;
		double *data;
	public:
		matrix(long n_rows, long n_cols, matrix_flags flags=uninitialized);
		matrix(long n_rows, long n_cols, double a11, ... );
		matrix(long n_rows, long n_cols, double *src);
		matrix(matrix const & m);
		~matrix( void );
		matrix & operator=(matrix const & m);
		matrix & operator+=(matrix const & m);
		matrix & operator-=(matrix const & m);
		matrix & operator*=(matrix const & m);
		matrix operator~(void);
		long n_rows(void);
		long n_cols(void);
		matrix inverse(void);
		void sub(matrix const & m, long row, long col);
		void sub(double s, long row, long col);
		matrix sub(long nrows, long ncols, long row, long col);
		void print(const char *name);
	};
	
	matrix operator+(matrix const & m);
	matrix operator-(matrix const & m);
	
	matrix operator+(matrix const & m1, matrix const & m2);
	matrix operator-(matrix const & m1, matrix const & m2);
	matrix operator*(matrix const & m1, matrix const & m2);
	matrix operator*(double const &s, matrix const & m);
	matrix operator*(matrix const & m, double const &s);
}

#endif
