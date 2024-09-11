#ifndef _TENSOR_H_
#define _TENSOR_H_

#include <iostream>
#include <vector>
#include "iterator.hpp"
class Vector;
class SparseVector;
class Matrix;
class SparseMatrix;

using vec2d = std::vector<std::vector<double>>;

Vector       operator + (const Vector &x,       const Vector &y);      
Vector       operator + (const Vector &x,       const SparseVector &y);
Vector       operator + (const SparseVector &x, const Vector &y);      
SparseVector operator + (const SparseVector &x, const SparseVector &y);
Matrix       operator + (const Matrix &x,       const Matrix &y);      
Matrix       operator + (const Matrix &x,       const SparseMatrix &y); 
Matrix       operator + (const SparseMatrix &x, const Matrix &y);  

Vector       operator - (const Vector &x,       const Vector &y);       
Vector       operator - (const Vector &x,       const SparseVector &y); 
Vector       operator - (const SparseVector &x, const Vector &y);       
SparseVector operator - (const SparseVector &x, const SparseVector &y); 
Matrix       operator - (const Matrix &x,       const Matrix &y);       
Matrix       operator - (const Matrix &x,       const SparseMatrix &y); 
Matrix       operator - (const SparseMatrix &x, const Matrix &y);   

double       operator * (const Vector &x,       const Vector &y);      
double       operator * (const SparseVector &x, const SparseVector &y);
double       operator * (const Vector &x,       const SparseVector &y);
double       operator * (const SparseVector &x, const Vector &y);      
Vector       operator * (double c,              const Vector &x);      
SparseVector operator * (double c,              const SparseVector &x);

Vector       operator * (const Matrix &mat,     const Vector &vec); 
SparseVector operator * (const Matrix &mat,     const SparseVector &vec);
SparseVector operator * (const SparseMatrix &x, const Vector &y); 
SparseVector operator * (const SparseMatrix &x, const SparseVector &y);
Matrix       operator * (const Matrix &x,       const Matrix &y);       
Matrix       operator * (double c,              const Matrix &x);       
SparseMatrix operator * (double c,              const SparseMatrix &x); 
SparseMatrix operator * (const SparseMatrix &x, const SparseMatrix &y); 
SparseMatrix operator * (const Matrix &x,       const SparseMatrix &y); 
SparseMatrix operator * (const SparseMatrix &x, const Matrix &y); 

Vector       operator / (const Vector &x,       double c); 
SparseVector operator / (const SparseVector &x, double c); 
Matrix       operator / (const Matrix &x,       double c); 
SparseMatrix operator / (const SparseMatrix &x, double c); 

std::ostream &operator << (std::ostream &stream, const Vector &vec);
std::ostream &operator << (std::ostream &stream, const SparseVector &vec);
std::ostream &operator << (std::ostream &stream, const Matrix &mat);
std::ostream &operator << (std::ostream &stream, const SparseMatrix &mat);

class Vector{ 

    protected:
        size_t n_size = 0;
        std::vector<double> values;
    public:
        Vector() = default;
        Vector(size_t n_size) : n_size(n_size) { }
        Vector(size_t n_size, double initial_value) : n_size(n_size) { values.resize(n_size, initial_value); }
        Vector(size_t n_size, std::vector<double> &values) : n_size(n_size), values(values) { }
 
        size_t getSize() const { return n_size; }
        double get(const int i) const { return values[i]; }
        const std::vector<double> &get() const { return values; }
        
        void operator += (const Vector &x);
        void operator += (const SparseVector &x);
        void operator -= (const Vector &x);
        void operator -= (const SparseVector &x);
        void operator *= (double c);
        void operator /= (double c);
        friend Vector operator + (const Vector &x, const Vector &y);
        friend Vector operator + (const Vector &x, const SparseVector &y);
        friend Vector operator + (const SparseVector &x, const Vector &y);
        friend Vector operator - (const Vector &x, const Vector &y);
        friend Vector operator - (const Vector &x, const SparseVector &y);
        friend Vector operator - (const SparseVector &x, const Vector &y);
        friend Vector operator * (double c, const Vector &x);
        friend double operator * (const Vector &x, const Vector &y);
        friend double operator * (const Vector &x, const SparseVector &y);
        friend double operator * (const SparseVector &x, const Vector &y);
        friend Vector operator * (const Matrix &mat, const Vector &vec);
        friend SparseVector operator * (const SparseMatrix &x, const Vector &y);
        friend Vector operator / (const Vector &x, double c);
        friend std::ostream &operator << (std::ostream &stream, const Vector &vec);
};

class SparseVector : public Vector { 
    
    private:
        std::vector<int> indices;
        using iterator = Iterator<std::vector<int>::iterator, std::vector<double>::iterator>;
        using const_iterator = Iterator<std::vector<int>::const_iterator, std::vector<double>::const_iterator>;
                
        iterator begin() { return iterator(indices.begin(), values.begin()); }
        iterator end() { return iterator(indices.end(), values.end()); }
        const_iterator begin() const { return const_iterator(indices.begin(), values.begin()); }
        const_iterator end() const { return const_iterator(indices.end(), values.end()); }

    public:
        SparseVector() = default;
        SparseVector(size_t n_size) : Vector(n_size) { }
        SparseVector(int n_size, std::vector<double> &values) : Vector(n_size, values) { }

        void push(int i, double value){
            indices.push_back(i);
            values.push_back(value);
        }

        friend Vector;
        //Denseベクトルを足してSparseになることが少ない想定で一部実装なし
        void operator += (const SparseVector &x);       
        void operator -= (const SparseVector &x);
        void operator *= (double c);
        void operator /= (double c);
        friend Vector       operator + (const Vector &x, const SparseVector &y);
        friend Vector       operator + (const SparseVector &x, const Vector &y);
        friend SparseVector operator + (const SparseVector &x, const SparseVector &y);
        friend SparseVector operator - (const SparseVector &x, const SparseVector &y);
        friend Vector       operator - (const Vector &x, const SparseVector &y);
        friend Vector       operator - (const SparseVector &x, const Vector &y);
        friend double       operator * (const SparseVector &x, const SparseVector &y);
        friend double       operator * (const Vector &x, const SparseVector &y);
        friend double       operator * (const SparseVector &x, const Vector &y);
        friend SparseVector operator * (double c, const SparseVector &x);
        friend SparseVector operator * (const Matrix &mat, const SparseVector &vec);
        friend SparseVector operator * (const SparseMatrix &x, const Vector &y);
        friend SparseVector operator * (const SparseMatrix &x, const SparseVector &y);
        friend SparseVector operator / (const SparseVector &x, double c);

        friend std::ostream &operator << (std::ostream &stream, const SparseVector &vec);
        void add_vec(const SparseVector &x, int sign);
};

class Matrix{ 
    protected: 
        size_t n_row, n_col;
        vec2d values;
    public:
        Matrix() = default;
        Matrix(size_t n_row, size_t n_col) : 
            n_row(n_row), n_col(n_col) { values.resize(n_row, std::vector<double>(n_col, 0)); }
        Matrix(size_t n_row, size_t n_col, vec2d values) : 
            n_row(n_row), n_col(n_col), values(values) { }

        int getRow() const { return n_row; }
        int getCol() const { return n_col; }
        double getValue(int r, int c) const { return values[r][c]; }

        void operator += (const Matrix &x);
        void operator += (const SparseMatrix &x);
        void operator -= (const Matrix &x);
        void operator -= (const SparseMatrix &x);
        //void operator *= (const Matrix &mat);
        void operator *= (double c);
        void operator /= (double c);
        friend Matrix       operator + (const Matrix &x, const Matrix &y);
        friend Matrix       operator + (const Matrix &x, const SparseMatrix &y);
        friend Matrix       operator + (const SparseMatrix &x, const Matrix &y);       
        friend Matrix       operator - (const Matrix &x, const Matrix &y);
        friend Matrix       operator - (const Matrix &x, const SparseMatrix &y);
        friend Matrix       operator - (const SparseMatrix &x, const Matrix &y);  
        friend Matrix       operator * (const Matrix &x, const Matrix &y);
        friend Matrix       operator * (double c, const Matrix &x);
        friend SparseMatrix operator * (const Matrix &x, const SparseMatrix &y);
        friend SparseMatrix operator * (const SparseMatrix &x, const Matrix &y);
        friend Vector       operator * (const Matrix &mat, const Vector &vec);
        friend SparseVector operator * (const Matrix &mat, const SparseVector &vec);
        friend Matrix       operator / (const Matrix &x, double c);
        friend std::ostream &operator << (std::ostream &stream, const Matrix &mat);
};

class SparseMatrix : public  Matrix {
private:
    std::vector<double> values;      // 非ゼロ要素の値
    std::vector<int> col_indices;    // 非ゼロ要素の列インデックス
    std::vector<int> row_pointers;   // 各行の非ゼロ要素の開始位置（行ポインタ）

public:
    // コンストラクタ
    SparseMatrix(int r, int c) : Matrix(r, c) {
        row_pointers.resize(r + 1, 0);  // 行数+1のサイズで初期化
    }

    void addValue(int row, int col, double value);
    void build() ;
    double getValue(int row, int col) const ;

    friend Matrix;
    void operator *= (double c);
    void operator /= (double c);
    friend Matrix       operator + (const Matrix &x, const SparseMatrix &y);
    friend Matrix       operator + (const SparseMatrix &x, const Matrix &y); 
    friend Matrix       operator - (const Matrix &x, const SparseMatrix &y);
    friend Matrix       operator - (const SparseMatrix &x, const Matrix &y);  
    friend SparseMatrix operator * (double c, const SparseMatrix &x);
    friend SparseMatrix operator * (const SparseMatrix &x, const SparseMatrix &y);
    friend SparseMatrix operator * (const Matrix &x, const SparseMatrix &y);
    friend SparseMatrix operator * (const SparseMatrix &x, const Matrix &y);
    friend SparseVector operator * (const SparseMatrix &x, const Vector &y);
    friend SparseVector operator * (const SparseMatrix &x, const SparseVector &y);
    friend SparseMatrix operator / (const SparseMatrix &x, double c);
    friend std::ostream &operator << (std::ostream &stream, const SparseMatrix &mat);

    //void operator += (const Matrix &x);
    //void operator += (const SparseMatrix &x);
    //void operator -= (const Matrix &x);
    //void operator -= (const SparseMatrix &x);
    //void operator *= (const Matrix &mat);
    //void operator *= (const SparseMatrix &mat);
    //friend Vector operator * (const Matrix &mat, const Vector &vec);
    //friend SparseVector operator * (const Matrix &mat, const SparseVector &vec);
    //friend SparseMatrix operator + (const SparseMatrix &x, const SparseMatrix &y); 
    //friend SparseMatrix operator - (const SparseMatrix &x, const SparseMatrix &y);     

};
#endif