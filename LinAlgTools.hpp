#ifndef LINALGTOOLS
#define LINALGTOOLS

#include <cmath>


template<typename ValueType, std::size_t dim>
class Vec;
template<typename ValueType, std::size_t m, std::size_t n>
class DenseMatrix;
template<typename ValueType, std::size_t m, std::size_t n>
class MatrixRow;
template<typename ValueType>
class MatrixView;


// Simple implementation of static sized vector class Vec that supports most standard
// mathematical operations like addition, subtraction, scalar multiplication, dot product
// 

template<typename ValueType, std::size_t dim>
class Vec {
    ValueType element[dim];
    public:
    using value_type = ValueType;
    using size_type = std::size_t;

    static constexpr std::size_t dimension = dim;

    Vec() {
        for (std::size_t i=0; i!=dim; ++i)
            element[i] = ValueType(0);
    }
    Vec(ValueType a) {
        for (std::size_t i=0; i!=dim; ++i)
            element[i] = a;
    }
    Vec(const std::initializer_list<ValueType> l) {
        //assert(l.size() == dim);
        std::copy(l.begin(), l.end(), element);
    }
    template<std::size_t m>
        Vec(MatrixRow<ValueType,m,dim> r) {
            ValueType *e = element;
            for (std::size_t i=0; i!=dim; ++i,++e)
                e = r[i];
        }
    template<typename _ValueType>
        Vec(Vec<_ValueType,dim> v) {
            for (std::size_t i=0; i!=dim; ++i)
                element[i] = v[i];
        }

    ValueType& operator[](std::size_t i) {
        return element[i];
    }
    const ValueType& operator[](std::size_t i) const {
        return element[i];
    }
    Vec<ValueType,dim>& operator+=(const Vec<ValueType,dim>& rhs) {
        for (std::size_t i=0; i!=dim; ++i)
            element[i] += rhs.element[i];
        return *this;
    }
    friend Vec<ValueType,dim> operator+(Vec<ValueType,dim> lhs,
            const Vec<ValueType,dim>& rhs) {
        lhs += rhs;
        return lhs;
    }
    Vec<ValueType,dim>& operator-=(const Vec<ValueType,dim>& rhs) {
        for (std::size_t i=0; i!=dim; ++i)
            element[i] -= rhs.element[i];
        return *this;
    }
    friend Vec<ValueType,dim> operator-(Vec<ValueType,dim> lhs,
            const Vec<ValueType,dim>& rhs) {
        lhs -= rhs;
        return lhs;
    }
    Vec<ValueType,dim>& operator*=(ValueType scal) {
        for (std::size_t i=0; i!=dim; ++i)
            element[i] *= scal;
        return *this;
    }
    friend Vec<ValueType,dim> operator*(ValueType scal, Vec<ValueType,dim> u) {
        return u*=scal;
    }
    friend Vec<ValueType,dim> operator*(Vec<ValueType,dim> u, ValueType scal) {
        return u*=scal;
    }
    Vec<ValueType,dim>& operator/=(ValueType scal) {
        for (std::size_t i=0; i!=dim; ++i)
            element[i] /= scal;
        return *this;
    }
    friend Vec<ValueType,dim> operator/(Vec<ValueType,dim> u, ValueType scal) {
        return u/=scal;
    }
    friend Vec<ValueType,dim> operator-(Vec<ValueType,dim> op) {
        op *= -1;
        return op;
    }
    friend ValueType prod(const Vec<ValueType,dim> &lhs,
            const Vec<ValueType,dim> &rhs) {
        ValueType product(0);
        for (std::size_t i=0; i!=dim; ++i)
            product += lhs.element[i]*rhs.element[i];
        return product;
    }
    friend ValueType norm_1(const Vec<ValueType,dim>& u) {
        using std::abs;
        ValueType norm(abs(u[0]));
        for (std::size_t i=1; i!=dim; ++i)
            norm += abs(u.element[i]);
        return norm;
    }
    friend ValueType norm_2(const Vec<ValueType,dim>& u) {
        using std::sqrt;
        ValueType product(u[0]*u[0]);
        for (std::size_t i=1; i!=dim; ++i)
            product += u.element[i]*u.element[i];
        return sqrt(product);
    }
    friend ValueType norm_max(const Vec<ValueType,dim>& u) {
        using std::abs;
        ValueType norm(abs(u[0]));
        for (std::size_t i=1; i!=dim; ++i)
            norm += abs(u.element[i]);
        return norm;
    }
    friend std::ostream& operator<<(std::ostream& stream,
            const Vec<ValueType,dim>& u) {
        ValueType const *b = u.element, *e = u.element + (dim-1);
        stream << '(';
        for (; b!=e; ++b)
            stream << *b << ',' << ' ';
        stream << *b << ')';
        return stream;
    }
    ValueType* begin() {
        return element;
    }
    ValueType* end() {
        return element+dim;
    }
    const ValueType* begin() const {
        return element;
    }
    const ValueType* end() const {
        return element+dim;
    }
    const ValueType* cbegin() const {
        return element;
    }
    const ValueType* cend() const {
        return element+dim;
    }
    constexpr std::size_t size() const {
        return dim;
    }
    friend void swap(Vec<ValueType,dim>& v1, Vec<ValueType,dim>& v2) {
        //std::cout << "Calling Vec::swap" << std::endl;
        using std::swap;
        for (std::size_t i=0; i!=dim; ++i) {
            swap(v1.element[i],v2.element[i]);
        }
    }
};

// Projects u onto line spanned by v
template<typename ValueType, std::size_t dim>
inline Vec<ValueType,dim> project(Vec<ValueType,dim> u,
        Vec<ValueType,dim> v) {
    return (prod(u,v)/norm_2(v))*v;
}

template<typename ValueType>
ValueType& first_element(ValueType& x) {return x;}
template<typename ValueType>
const ValueType& first_element(const ValueType& x) {return x;}
template<typename ValueType, std::size_t dim>
ValueType& first_element(Vec<ValueType,dim>& x) {return x[0];}
template<typename ValueType, std::size_t dim>
const ValueType& first_element(const Vec<ValueType,dim>& x) {return x[0];}

// Column-major for now
template<typename ValueType, std::size_t m, std::size_t n>
class DenseMatrix {
    private:
        ValueType element[m*n];
    public:
        using value_type = ValueType;
        using size_type = std::size_t;

        DenseMatrix() {
            for (std::size_t i=0; i!=m*n; ++i)
                element[i] = 0;
        }
        DenseMatrix(ValueType a) {
            for (std::size_t i=0; i!=m*n; ++i)
                element[i] = a;
        }
        DenseMatrix(const std::initializer_list<ValueType> l) {
            //assert(l.size() == m*n);
            std::copy(l.begin(), l.end(), element);
        }

        MatrixRow<ValueType,m,n> operator[](std::size_t i) {
            return MatrixRow<ValueType,m,n>(*this,i);
        }
        const MatrixRow<ValueType,m,n> operator[](std::size_t i) const {
            return MatrixRow<ValueType,m,n>(*this,i);
        }

        Vec<ValueType,m> column(std::size_t i) {
            Vec<ValueType,m> u;
            std::copy(element+i*m,element+(i+1)*m,u.begin());
            return u;
        }
        void column(Vec<ValueType,m> u, std::size_t i) {
            std::copy(u.begin(),u.end(),element+i*m);
        }

        ValueType* data() {
            return element;
        }

        friend std::ostream& operator<<(std::ostream& stream,
                const DenseMatrix<ValueType,m,n>& A) {
            ValueType const *i, *b = A.element, *e = b + m*(n-1), *final = b + m; 
            for (;b != final; ++b,++e) {
                stream << '[';
                for (i=b; i!=e; i+=m)
                    stream << *i << ',' << ' ';
                stream << *i << ']' << std::endl;
            }
            return stream;
        }

        template<std::size_t k>
            friend DenseMatrix<ValueType,m,n> operator*(DenseMatrix<ValueType,m,k>& A,DenseMatrix<ValueType,k,n>& B) {
                DenseMatrix<ValueType,m,n> C(0);
                for (size_type i=0; i!=m; ++i)
                    for (size_type j=0; j!=n; ++j)
                        for (size_type l=0; l!=k; ++l)
                            C[i][j] += A[i][l]*B[l][j];
                return C;
            }
};

template<typename ValueType, std::size_t m, std::size_t n, std::size_t k>
void column_copy(DenseMatrix<ValueType,m,n>& A, std::size_t nstart, std::size_t nend,DenseMatrix<ValueType,m,k>& B, std::size_t kstart) {
    for (std::size_t j=0; j!=nend-nstart; ++j) {
        for (std::size_t i=0; i!=m; ++i) {
            B[i][kstart+j] = A[i][nstart+j];
            // std::cout << "i=" << i << ", j=" << j << std::endl;
        }
    }
}

template<typename ValueType, std::size_t m, std::size_t n, std::size_t k>
void column_copy_rev(DenseMatrix<ValueType,m,n>& A, std::size_t nstart, std::size_t nend,DenseMatrix<ValueType,m,k>& B, std::size_t kstart) {
    for (std::size_t j=0; j!=nend-nstart; ++j) {
        for (std::size_t i=0; i!=m; ++i) {
            B[i][kstart-j] = A[i][nstart+j];
            // std::cout << "i=" << i << ", j=" << j << std::endl;
        }
    }
}

template<typename ValueType, std::size_t m, std::size_t n, std::size_t k, std::size_t l>
void column_fuse(DenseMatrix<ValueType,m,n>& A, std::size_t nstart, std::size_t nend,DenseMatrix<ValueType,m,k>& B, std::size_t kstart, std::size_t kend, DenseMatrix<ValueType,m,l>& C) {
    for (std::size_t j=0; j!=nend-nstart; ++j) {
        for (std::size_t i=0; i!=m; ++i) {
            C[i][j] = A[i][j+nstart];
            // std::cout << "i=" << i << ", j=" << j << std::endl;
        }
    }
    for (std::size_t j=0; j!=kend-kstart; ++j) {
        for (std::size_t i=0; i!=m; ++i) {
            C[i][j+nend-nstart] = B[i][j+kstart];
            // std::cout << "i=" << i << ", j=" << j << std::endl;
        }
    }
}


template<typename ValueType, std::size_t m, std::size_t n>
class MatrixRow {
    private:
        ValueType* matptr;
        std::size_t rownr;
    public:
        MatrixRow() : matptr(nullptr) {}
        MatrixRow(DenseMatrix<ValueType,m,n>& mat, std::size_t _rownr) : matptr(mat.data()), rownr(_rownr) {};

        ValueType& operator[](std::size_t i) {
            return matptr[rownr+i*m];
        }
        const ValueType& operator[](std::size_t i) const {
            return matptr[rownr+i*m];
        }
};

template<typename ValueType>
class MatrixView{
    private:
        ValueType *ptr;
        std::size_t Cols,Col_stride,Rows,Row_stride;
    public:
        template<std::size_t m,std::size_t n>
            MatrixView(DenseMatrix<ValueType,m,n> A) {
                ptr = A.data();
                Cols = m; Rows = n;
                Col_stride = 1; Row_stride = m;
            }
        template<std::size_t m,std::size_t n>
            MatrixView(DenseMatrix<ValueType,m,n> A, std::size_t mstart,
                    std::size_t mend, std::size_t nstart, std::size_t nend) {
                ptr = A.data()+nstart*m+mstart;
                Cols = mend-mstart+1; Rows = nend-nstart+1;
                Col_stride = 1; Row_stride = m;
            }
        template<std::size_t m,std::size_t n>
            MatrixView(DenseMatrix<ValueType,m,n> A, std::size_t mstart,
                    std::size_t mstep, std::size_t mend, std::size_t nstart,
                    std::size_t nstep, std::size_t nend) {
                ptr = A.data()+nstart*m+mstart;
                Cols = mend-mstart+1; Rows = nend-nstart+1;
                Col_stride = mstep; Row_stride = m*nstep;
            }
};


////////////////////////////////////
/////////// ALGORITHMS /////////////
////////////////////////////////////


//  QR-factorization algorithms : A = QR ( A-mxn, Q-mxn, R-nxn )
template<typename ValueType, std::size_t m, std::size_t n>
void QR_Gram_Schmidt(DenseMatrix<ValueType,m,n>& A,
        DenseMatrix<ValueType,m,n>& Q,
        DenseMatrix<ValueType,n,n>& R) {
    // assert(m >= n)
    for (std::size_t i=0; i!=n; ++i) {
        for (std::size_t j=0; j!=i; ++j) {
            //R[j,i] = Q[:,j]*A[:,i]
            R[j][i] = 0;
            for (std::size_t k=0; k!=m; ++k)
                R[j][i] += Q[k][j]*A[k][i];
        }
        //Q[:,i] = A[:,i] - SUM_{j=0,...,i-1} R[j,i]Q[:,j]
        for (std::size_t k=0; k!=m; ++k) {
            Q[k][i] = A[k][i];
            for (std::size_t j=0; j!=i; ++j) {
                Q[k][i] -= R[j][i]*Q[k][j];
            } 
        }
        //R[i,i] = ||Q[:,i]||
        R[i][i] = 0;
        for (std::size_t k=0; k!=m; ++k)
            R[i][i] += Q[k][i]*Q[k][i];
        using std::sqrt;
        R[i][i] = sqrt(R[i][i]);
        //Q[:,i] = Q[:,i]/R[i,i]
        for (std::size_t k=0; k!=m; ++k)
            Q[k][i] = Q[k][i]/R[i][i];
    }
}

// Function for reducing AX = B to A'X=B' where A' is upper triangular (without solving for X)
// A is mxn-matrix and B is mxk (so X is nxk)
template<typename ValueType,std::size_t m,std::size_t n,std::size_t k>
void triangular(DenseMatrix<ValueType,m,n>& A,DenseMatrix<ValueType,m,k>& B) {
    ValueType a(0);
    for (std::size_t i=0; i!=n; ++i) {
        for (std::size_t j=i+1; j!=m; ++j) {
            a = A[j][i]/A[i][i];
            A[j][i] = 0;
            for (std::size_t l=i+1; l!=n; ++l) {
                A[j][l] -= A[i][l]*a;
            }
            for (std::size_t l=0; l!=k; ++l) {
                B[j][l] -= B[i][l]*a;
            }
        }
    }
}

// Perform backwards substitution on AX=B where A is upper triangular
// A is mxm, B is mxk. If X is mxk the full substitution is done but if X is nxk, n<m, only partial substitution is done (starting from the back)
template<typename ValueType,std::size_t m,std::size_t n,std::size_t k>
void back_sub(DenseMatrix<ValueType,m,m>& A,DenseMatrix<ValueType,m,k>& B,DenseMatrix<ValueType,n,k>& X, std::size_t kx) {
    ValueType temp(0);
    for (std::size_t i=1; i!=n+1; ++i) {
        // std::cout << "i=" << i << std::endl;
        for (std::size_t l=0; l!=kx; ++l) {
            // std::cout << " l=" << l << std::endl;
            temp = B[m-i][l];
            // std::cout << " temp=" << temp << std::endl;
            for (std::size_t j=1; j!=i; ++j) {
                // std::cout << "  j=" << j << std::endl;
                temp -= A[m-i][m-i+j]*X[m-i+j][l];
                // std::cout << "  temp=" << temp << std::endl;
            }
            X[m-i][l] = temp/A[m-i][m-i];
        }
    }
}

template<typename ValueType, std::size_t dim>
ValueType max(Vec<ValueType,dim> v) {
    ValueType m = v[0];
    for (auto i = v.begin()+1; i != v.end(); ++i)
        if ( m < *i)
            m = *i;
    return m;
}

template<typename ValueType, std::size_t dim>
ValueType min(Vec<ValueType,dim> v) {
    ValueType m = v[0];
    for (auto i = v.begin()+1; i != v.end(); ++i)
        if ( m > *i)
            m = *i;
    return m;
}

// Pass by value and remove B temp?
template<typename ValueType, std::size_t m>
ValueType det(DenseMatrix<ValueType,m,m>& A) {
    DenseMatrix<ValueType,m,0> B;
    triangular(A,B);
    ValueType ans = 1;
    for (std::size_t i=0; i!=m; ++i)
        ans *= A[i][i];
    return ans;
}

template<typename ValueType, std::size_t m, std::size_t n>
ValueType Gram_det(DenseMatrix<ValueType,m,n>& A) {
    DenseMatrix<ValueType,n,n>& G(0);
    for (std::size_t i=0; i!=n; ++i)
        for (std::size_t j=0; j!=n; ++j)
            for (std::size_t k=0; k!=m; ++k)
                G[i][j] += A[i][k]*A[j][k];
    return det(G);
}

// template<typename ValueType, std::size_t m, std::size_t n, std::size_t k>
// void flag_reduce(DenseMatrix<ValueType,m,n>& V, DenseMatrix<ValueType,m,n>& W,DenseMatrix<ValueType,m,n>& U, Vec<std::size_t,k> lambdadim) {
//   DenseMatrix<ValueType,m,n> A, B, X;
//   std::size_t dim = 0;
//   for (std::size_t l=0; l!=k; ++l) {
//     column_copy(V,dim,dim+lambdadim[l],B,0);
//     column_copy(W,0,n-dim-1,A,0);
//     column_copy_rev(V,0,dim,A,n-1);
//     triangle(A,B);
//     back_sub(A,B,X,lambdadim[l]);
//     dim += lambdadim[l];

//     column_copy();
//   }
// }


#endif
