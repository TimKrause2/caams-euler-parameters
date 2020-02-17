#include <cstdlib>
#include <cstdarg>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "matrix.hpp"

namespace caams
{
	matrix::matrix(long n_rows, long n_cols, matrix_flags flags) {
		rows = n_rows;
		cols = n_cols;
		long N=n_rows*n_cols;
		data = new double[N];
		if(flags&init){
			double *d = data;
			if(flags&identity){
				if(rows!=cols){
					std::cout << "matrix::matrix : identity requested but rows not equal columns" << std::endl;
					std::exit(1);
				}
				for(long r=0;r<rows;r++){
					for(long c=0;c<cols;c++){
						if(r==c){
							*(d++) = 1.0;
						}else{
							*(d++) = 0.0;
						}
					}
				}
			}else{
				for(long i=0;i<N;i++){
					*(d++)=0.0;
				}
			}
		}
	}
	
	matrix::matrix(long n_rows, long n_cols, double a11, ...){
		rows = n_rows;
		cols = n_cols;
		long N=n_rows*n_cols;
		data = new double[N];
		va_list args;
		va_start(args, a11);
		double *d = data;
		*(d++)=a11;
		for(long i=1;i<N;i++){
			*(d++) = va_arg(args, double);
		}
		va_end(args);
	}
	
	matrix::matrix(long n_rows, long n_cols, double *src){
		rows = n_rows;
		cols = n_cols;
		long N=rows*cols;
		data = new double[N];
		double *d = data;
		for(long i=N;i;i--){
			*(d++) = *(src++);
		}
	}
	
	matrix::matrix(matrix const & m){
		rows = m.rows;
		cols = m.cols;
		long N=rows*cols;
		data = new double[N];
		double *d = data;
		double *s = m.data;
		for(long i=0;i<N;i++){
			*(d++) = *(s++);
		}
	}
	
	matrix::~matrix(void){
		delete [] data;
	}
	
	matrix & matrix::operator=(matrix const & m){
		if(rows!=m.rows || cols!=m.cols){
			std::cout << "matrix::operator= : unequal matrix dimensions" << std::endl;
			std::exit(1);
		}
		double *d = data;
		double *s = m.data;
		long N=rows*cols;
		for(long i=0;i<N;i++){
			*(d++)=*(s++);
		}
		return *this;
	}
	
	matrix & matrix::operator+=(matrix const & m){
		if(rows!=m.rows || cols!=m.cols){
			std::cout << "matrix::operator+= : unequal matrix dimensions" << std::endl;
			std::exit(1);
		}
		double *d = data;
		double *s = m.data;
		long N=rows*cols;
		for(long i=0;i<N;i++){
			*(d++)+=*(s++);
		}
		return *this;
	}
	
	matrix & matrix::operator-=(matrix const & m){
		if(rows!=m.rows || cols!=m.cols){
			std::cout << "matrix::operator-= : unequal matrix dimensions" << std::endl;
			std::exit(1);
		}
		double *d = data;
		double *s = m.data;
		long N=rows*cols;
		for(long i=0;i<N;i++){
			*(d++)-=*(s++);
		}
		return *this;
	}
	
	matrix & matrix::operator*=(matrix const &m){
		if(rows!=cols || rows!=m.rows || cols!=m.cols){
			std::cout << "matrix::operator*= : inconsistent matrix dimensions" << std::endl;
			std::exit(1);
		}
		*this = *this * m;
		return *this;
	}
	
	matrix matrix::operator~(void){
		matrix res(cols,rows);
		double *d = res.data;
		double *s = data;
		for(long r=0;r<rows;r++){
			for(long c=0;c<cols;c++){
				*d = *(s++);
				d += res.cols;
			}
			d -= res.cols*res.rows-1;
		}
		return res;
	}
	
	long matrix::n_rows(void){
		return rows;
	}
	
	long matrix::n_cols(void){
		return cols;
	}
	
	matrix matrix::inverse(void){
		if(rows!=cols){
			std::cout << "matrix::inverse : rows not equal columns" << std::endl;
			std::exit(1);
		}
		// create identity matrix
		//dmatrixIdentity( inv );
		matrix inv(rows,cols,init_identity);

		// make a copy of the original
		//dmatrix *m_copy = dmatrixNew( m->n_rows, m->n_cols );
		//dmatrixCopy( m, m_copy );
		matrix m_copy( *this );
		
		// perform row operations to obtain echelon form
		long pivot;
		double *pivot0 = m_copy.data;
		double *row_inv0 = inv.data;
		for(pivot=0;pivot<(rows-1);pivot++){
			// find the row with the largest absolute value
			double *col_data = pivot0;
			double *col_inv_data = row_inv0;
			long row;
			double d_max = 0.0;
			long row_max;
			double *row_max_data;
			double *row_inv_max_data;
			for(row=pivot;row<rows;row++){
				double d_abs = std::fabs(*col_data);
				if( d_abs > d_max ){
					d_max = d_abs;
					row_max = row;
					row_max_data = col_data;
					row_inv_max_data = col_inv_data;
				}
				col_data += cols;
				col_inv_data += inv.cols;
			}
			// use the maximum row to reduce the other rows under the pivot
			col_data = pivot0;
			col_inv_data = row_inv0;
			for(row=pivot;row<rows;row++){
				if(row!=row_max){
					long col;
					double *data = col_data;
					double *data_reduce = row_max_data;
					double alpha = data[0]/data_reduce[0];
					for(col=pivot;col<cols;col++){
						*data -= alpha * *data_reduce;
						data++;
						data_reduce++;
					}
					col_data[0] = 0.0;
					data = col_inv_data;
					data_reduce = row_inv_max_data;
					for(col=0;col<inv.cols;col++){
						*data -= alpha * *data_reduce;
						data++;
						data_reduce++;
					}
				}
				col_data += cols;
				col_inv_data += inv.cols;
			}

			// exchange the maximum row with the pivot
			if( row_max!=pivot ){
				double *src0 = pivot0;
				double *src1 = row_max_data;
				long col;
				double d_temp;
				for(col=pivot;col<cols;col++){
					d_temp = *src0;
					*src0 = *src1;
					*src1 = d_temp;
					src0++;
					src1++;
				}
				src0 = row_inv0;
				src1 = row_inv_max_data;
				for(col=0;col<inv.cols;col++){
					d_temp = *src0;
					*src0 = *src1;
					*src1 = d_temp;
					src0++;
					src1++;
				}
			}
			pivot0 += cols + 1;
			row_inv0 += inv.cols;
		}
		
		//printf("Matricis after row reduction:\n");
		//dmatrixPrint( m_copy, "m_copy" );
		//dmatrixPrint( inv, "inv" );
		//std::cout << "Matricis after row reduction:" << std::endl;
		//m_copy.print("m_copy");
		//inv.print("inv");
		
		// normalize the diagonal elements
		pivot0 = m_copy.data;
		row_inv0 = inv.data;
		for(pivot=0;pivot<rows;pivot++){
			// normalize this row
			double alpha = 1.0 / *pivot0;
			double *data = pivot0;
			long col;
			for(col=pivot;col<cols;col++){
				*(data++) *= alpha;
			}
			*pivot0 = 1.0;
			data = row_inv0;
			for(col=0;col<inv.cols;col++){
				*(data++) *= alpha;
			}
			pivot0 += cols + 1;
			row_inv0 += inv.cols;
		}

		//printf("Matricis after diagonal normalization:\n");
		//dmatrixPrint( m_copy, "m_copy" );
		//dmatrixPrint( inv, "inv" );
		//std::cout << "Matricis after diagonal normalization:" << std::endl;
		//m_copy.print("m_copy");
		//inv.print("inv");
		
		// perform row operations to obtain identity
		pivot0 = m_copy.data + m_copy.cols*m_copy.rows - 1;
		row_inv0 = inv.data + inv.cols*(inv.rows-1);
		for(pivot=cols;pivot>1;pivot--){
			double *col_data = pivot0 - cols;
			double *col_inv_data = row_inv0 - inv.cols;
			long row;
			for(row=pivot-1;row;row--){
				double alpha = *col_data;
				// reduce the row in the original
				*col_data -= alpha; // * *pivot0;
				// reduce the row in the inverse
				double *data = col_inv_data;
				double *data_reduce = row_inv0;
				long col;
				for(col=0;col<inv.cols;col++){
					*(data++) -= *(data_reduce++) * alpha;
				}
				col_data -= cols;
				col_inv_data -= inv.cols;
			}
			pivot0 -= cols + 1;
			row_inv0 -= inv.cols;
		}

		//printf("Matricis after upper row reduction:\n");
		//dmatrixPrint( m_copy, "m_copy" );
		//dmatrixPrint( inv, "inv" );
		//std::cout << "Matricis after upper row reduction:" << std::endl;
		//m_copy.print("m_copy");
		//inv.print("inv");

		// free the copy
		//dmatrixFree( m_copy );
		
		return inv;
	}
	
	void matrix::sub(matrix const & m, long row, long col){
		row--;
		col--;
		if( (m.rows+row)>rows || (m.cols+col)>cols ){
			std::cout << "matrix::sub : sub matrix can not fit" << std::endl;
			std::exit(1);
		}
		double *d = data + row*cols + col;
		double *s = m.data;
		for(long r=0;r<m.rows;r++){
			for(long c=0;c<m.cols;c++){
				*(d++) = *(s++);
			}
			d+=cols-m.cols;
		}
	}
	
	void matrix::sub(double s, long row, long col){
		row--;
		col--;
		if(row>=rows || col>=cols){
			std::cout << "matrix::sub : row or col outside of matrix" << std::endl;
			std::exit(1);
		}
		data[row*cols+col] = s;
	}
	
	matrix matrix::sub(long nrows, long ncols, long row, long col){
		row--;
		col--;
		if( (row+nrows)>rows || (col+ncols)>cols ){
			std::cout << "matrix::sub : sub matrix beyond matrix size" << std::endl;
			std::exit(1);
		}
		matrix res(nrows,ncols);
		double *d = res.data;
		double *s = &data[row*cols+col];
		for(long r=0;r<nrows;r++){
			for(long c=0;c<ncols;c++){
				*(d++) = *(s++);
			}
			s+=cols-ncols;
		}
		return res;
	}
	
	void matrix::print(const char *m){
		std::cout << m << std::endl;
		double *d = data;
		for(long r=0;r<rows;r++){
			for(long c=0;c<cols;c++){
                std::cout << std::setprecision(15) << *(d++);
				if(c<(cols-1)){
					std::cout << " ";
				}
			}
			std::cout << std::endl;
		}
	}
	
	matrix operator+(matrix const & m){
		return m;
	}
	
	matrix operator-(matrix const & m){
		return -1.0*m;
	}
	
	matrix operator+(matrix const & m1, matrix const & m2){
		if(m1.rows!=m2.rows || m1.cols!=m2.cols){
			std::cout << "operator+ : matricis are different dimensions" << std::endl;
			std::exit(1);
		}
		matrix r(m1.rows,m1.cols);
		double *d = r.data;
		double *s1 = m1.data;
		double *s2 = m2.data;
		for(long i=m1.rows*m1.cols;i;i--){
			*(d++) = *(s1++) + *(s2++);
		}
		return r;
	}
	
	matrix operator-(matrix const & m1, matrix const & m2){
		if(m1.rows!=m2.rows || m1.cols!=m2.cols){
			std::cout << "operator- : matricis are different dimensions" << std::endl;
			std::exit(1);
		}
		matrix r(m1.rows,m1.cols);
		double *d = r.data;
		double *s1 = m1.data;
		double *s2 = m2.data;
		for(long i=m1.rows*m1.cols;i;i--){
			*(d++) = *(s1++) - *(s2++);
		}
		return r;
	}
	
	matrix operator*(matrix const &m1, matrix const &m2){
		if(m1.cols!=m2.rows){
			std::cout << "operator*(matrix,matrix) : inconsistent matrix dimensions" << std::endl;
			std::exit(1);
		}
		matrix res(m1.rows,m2.cols);
		double *d = res.data;
		double *s1 = m1.data;
		double *s2 = m2.data;
		for(long r=0;r<res.rows;r++){
			for(long c=0;c<res.cols;c++){
				*d = 0.0;
				for(long i=0;i<m1.cols;i++){
					*d += *s1 * *s2;
					s1++;
					s2+=m2.cols;
				}
				d++;
				s1-=m1.cols;
				s2-=m2.cols*m2.rows-1;
			}
			s1+=m1.cols;
			s2=m2.data;
		}
		return res;
	}
	
	matrix operator*(double const &t, matrix const & m){
		matrix res(m.rows,m.cols);
		double *d = res.data;
		double *s = m.data;
		for(long i=m.rows*m.cols;i;i--){
			*(d++) = *(s++) * t;
		}
		return res;
	}
	
	matrix operator*(matrix const &m, double const &t){
		return t*m;
	}
	
}
