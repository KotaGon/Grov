#include "tensor.hpp"
//=======================
// Vector 
//=======================
void Vector::operator+=(const Vector &x){ 
    assert(n_size == x.n_size);
    for (size_t i = 0; i < n_size; i++)
        values[i] += x.values[i];    
}

void Vector::operator+=(const SparseVector &x){ 
    for(const auto [indice, value] : x) 
        values[indice] += value; 
}

void Vector::operator-=(const Vector &x){
    assert(n_size == x.n_size);
    for(size_t i = 0; i < n_size; i++)  
        values[i] -= x.values[i];
}

void Vector::operator -= (const SparseVector &x){ 
       for(const auto [indice, value] : x) 
        values[indice] -= value; 
}

void Vector::operator*=(double c){
    for(auto &value : values)
        value *= c;   
}

void Vector::operator/=(double c){
    *this *= 1.0 / c;
}

void SparseVector::operator+=(const SparseVector &x){ 
    add_vec(x, 1);
}

void SparseVector::operator -=(const SparseVector &x){ 
    add_vec(x, -1);
}

void SparseVector::operator *= (double c){ 
    for(auto &value : values)
        value *= c;
}

void SparseVector::operator /= (double c){ 
    *this *= 1.0 / c;
}

Vector operator + (const Vector &x, const Vector &y){ 
    Vector z = x;
    z += y;
    return z;
}

Vector operator + (const Vector &x, const SparseVector &y){ 
    Vector z = x;
    z += y;
    return z;
}

Vector operator + (const SparseVector &x, const Vector &y){ 
    Vector z = y;
    z += x;
    return z;
}

SparseVector operator + (const SparseVector &x, const SparseVector &y){ 
    SparseVector z = x;
    z += y;
    return z;
}

Vector operator - (const Vector &x, const Vector &y){ 
    return x + (-1) * y;
}

Vector operator - (const Vector &x, const SparseVector &y){
    Vector z = x;
    z -= y;
    return z;
}

Vector operator - (const SparseVector &x, const Vector &y){ 
    Vector z = (-1) * y;
    z += x;
    return z;
}

SparseVector operator - (const SparseVector &x, const SparseVector &y){ 
    SparseVector z = x;
    z -= y;
    return z;
}

double operator * (const Vector &x, const Vector &y){ 
    assert(x.n_size == y.n_size);
    double dot = 0;
    for (size_t i = 0; i < x.n_size; i++)
        dot += x.values[i] * y.values[i];
    return dot;    
}

double operator * (const SparseVector &x, const SparseVector &y){ 
    assert(x.n_size == y.n_size);
    double dot = 0;
    size_t i = 0, j = 0;
    while (i < y.indices.size() && j < x.indices.size()) {
        if (y.indices[i] == x.indices[j]) {
            // 同じインデックスがあれば、値を足す
            dot += y.values[i] * x.values[j];
            ++i; ++j;
        } else if (y.indices[i] < x.indices[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    return dot;  
}

double operator * (const Vector &x, const SparseVector &y){ 
    double dot = 0;
    for(const auto [indice, value] : y)
        dot += x.values[indice] * value;
    return dot;    
}

double operator * (const SparseVector &x, const Vector &y) { 
    return y * x;
}

Vector operator * (double c, const Vector &x){ 
    Vector y = x;
    y *= c;
    return y;
}

SparseVector operator * (double c, const SparseVector &x){ 
    SparseVector y = x;
    y *= c;
    return y;
}

Vector operator / (const Vector &x, double c){ 
    Vector y = x;
    y /= c;
    return y;
}

SparseVector operator / (const SparseVector &x, double c){ 
    SparseVector y = x;
    y /= c;
    return y;
}

std::ostream &operator << (std::ostream &stream, const Vector &vec){
    for (size_t i = 0; i < vec.getSize(); i++)
        stream << "Vec" << i << ": " << vec.get(i) << std::endl;
    return stream;
}

std::ostream &operator << (std::ostream &stream, const SparseVector &vec){  
    for(const auto [indice, value] : vec)
        stream << "Vec" << indice << ": " << value << std::endl;
    return stream;
}

void SparseVector::add_vec(const SparseVector &x, int sign){ 
    assert(n_size == x.n_size);  // サイズが同じであることを確認

    std::vector<double> result_values;
    std::vector<int> result_indices;

    size_t i = 0, j = 0;

    // マージ処理
    while (i < indices.size() && j < x.indices.size()) {
        if (indices[i] == x.indices[j]) {
            // 同じインデックスがあれば、値を足す
            double sum = values[i] + sign * x.values[j];
            if (sum != 0.0) {
                result_values.push_back(sum);
                result_indices.push_back(indices[i]);
            }
            ++i;
            ++j;
        } else if (indices[i] < x.indices[j]) {
            // thisのインデックスが小さい場合
            result_values.push_back(values[i]);
            result_indices.push_back(indices[i]);
            ++i;
        } else {
            // xのインデックスが小さい場合
            result_values.push_back(sign * x.values[j]);
            result_indices.push_back(x.indices[j]);
            ++j;
        }
    }

    // 残りの要素を追加（this側）
    while (i < indices.size()) {
        result_values.push_back(values[i]);
        result_indices.push_back(indices[i]);
        ++i;
    }

    // 残りの要素を追加（x側）
    while (j < x.indices.size()) {
        result_values.push_back(sign * x.values[j]);
        result_indices.push_back(x.indices[j]);
        ++j;
    }

    // 結果を現在のベクトルに反映
    values = result_values;
    indices = result_indices;   
}


//=======================
// Matrix 
//=======================
void Matrix::operator += (const Matrix &x){ 
    assert(n_row == x.n_row && n_col == x.n_col);
    for (size_t i = 0; i < n_row; i++)
        for (size_t j = 0; j < n_col; j++)
            values[i][j] += x.values[i][j];
}

void Matrix::operator += (const SparseMatrix &x){ 
    for (size_t i = 0; i < n_row; i++)
        for (size_t j = x.row_pointers[i]; j < x.row_pointers[i + 1]; ++j)
            values[i][x.col_indices[j]] += x.values[j];
}

void Matrix::operator -= (const Matrix &x){ 
    assert(n_row == x.n_row && n_col == x.n_col);
    for (size_t i = 0; i < n_row; i++)
        for (size_t j = 0; j < n_col; j++)
            values[i][j] -= x.values[i][j];
}

void Matrix::operator -= (const SparseMatrix &x){ 
    for (size_t i = 0; i < n_row; i++)
        for (size_t j = x.row_pointers[i]; j < x.row_pointers[i + 1]; ++j)
            values[i][x.col_indices[j]] -= x.values[j];
}

void Matrix::operator *= (double c) { 
    for (size_t i = 0; i < n_row; i++)
        for (size_t j = 0; j < n_col; j++)
            values[i][j] *= c; 
}

void SparseMatrix::operator *= (double c) { 
    for(auto &value : values) 
        value *= c;
}

void Matrix::operator /= (double c) { 
    *this *= 1.0 / c;
}

void SparseMatrix::operator /= (double c) { 
    *this *= 1.0 / c;
}

Matrix operator + (const Matrix &x, const Matrix &y){ 
    assert(x.n_row == y.n_row && x.n_col == y.n_col);
    Matrix z = x;
    z += y;
    return z;
}

Matrix operator + (const Matrix &x, const SparseMatrix &y){ 
    Matrix z = x;
    z += y;
    return z;
}

Matrix operator + (const SparseMatrix &x, const Matrix &y){
    Matrix z = y;
    z += x;
    return z;
}

Matrix operator - (const Matrix &x, const Matrix &y){ 
    assert(x.n_row == y.n_row && x.n_col == y.n_col);
    Matrix z = x;
    z -= y;
    return z;
}

Matrix operator - (const Matrix &x, const SparseMatrix &y){ 
    Matrix z = x;
    z -= y;
    return z;
}

Matrix operator - (const SparseMatrix &x, const Matrix &y){ 
    Matrix z = (-1) * y;
    z += x;
    return z;
}

Matrix operator * (const Matrix &x, const Matrix &y){ 
    assert(x.n_col == y.n_row);
    Matrix z = Matrix(x.n_row, y.n_col);
    for (size_t i = 0; i < x.n_row; i++)
        for (size_t j = 0; j < y.n_col; j++)
            for (size_t k = 0; k < x.n_col; k++)
                z.values[i][j] += x.values[i][k] * y.values[k][j];
    return z;
}

Matrix operator * (double c, const Matrix &x){
    Matrix z = x;
    z *= c;
    return z;
}

SparseVector operator * (const SparseMatrix &mat, const Vector &vec){ 
    assert(mat.n_col == vec.n_size);
    SparseVector y;
    for (size_t i = 0; i < mat.n_row; i++)
    {
        double element = 0.0;
        for (size_t j = mat.row_pointers[i]; j < mat.row_pointers[i + 1]; ++j)
            element += mat.values[j] * vec.values[mat.col_indices[j]];
        if(element != 0.0) 
            y.push(i, element);
    }
    return y;
}

SparseVector operator * (const SparseMatrix &mat, const SparseVector &vec){ 
    assert(mat.n_col == vec.n_size);
    SparseVector y;
    std::vector<double> vec_values(mat.n_col, 0.0);

    for(const auto [indice, value] : vec)   
        vec_values[indice] = value;

    for (size_t i = 0; i < mat.n_row; i++)
    {
        double element = 0.0;
        for (size_t j = mat.row_pointers[i]; j < mat.row_pointers[i + 1]; ++j)
            element += mat.values[j] * vec_values[mat.col_indices[j]];
        if(element != 0.0) 
            y.push(i, element);        
    }
    
    return y;
}

SparseMatrix operator * (const SparseMatrix &sparse, const Matrix &dense){ 
    assert(sparse.n_col == dense.n_row);  // 行列のサイズが適切であることを確認

    int rows = sparse.n_row;
    int cols = dense.n_col;
    SparseMatrix result(rows, cols);  // 結果の疎行列

    // 各行についてループ
    for (int i = 0; i < rows; ++i) {
        // 疎行列のi行目の非ゼロ要素に対してループ
        for (int j = sparse.row_pointers[i]; j < sparse.row_pointers[i + 1]; ++j) {
            int a_col = sparse.col_indices[j];
            double a_value = sparse.values[j];

            // 密行列の各列に対して計算
            for (int k = 0; k < cols; ++k) {
                double result_value = a_value * dense.values[a_col][k];
                if (result_value != 0.0) {
                    result.addValue(i, k, result_value);
                }
            }
        }
    }

    result.build();  // 行ポインタを構築
    return result;
}

SparseMatrix operator*(const Matrix& dense, const SparseMatrix& sparse) {
    assert(dense.n_col == sparse.n_row);

    int rows = dense.n_row;
    int cols = sparse.n_col;
    SparseMatrix result(rows, cols);

    // 一時的な結果を保持するハッシュマップ（列インデックスをキー、値が結果）
    std::unordered_map<int, double> temp;

    // 密行列の各行についてループ
    for (int i = 0; i < rows; ++i) {
        temp.clear();  // ハッシュマップをクリア

        // 疎行列の各行に対してループ
        for (int j = 0; j < sparse.n_row; ++j) {
            // 疎行列のj行目の非ゼロ要素に対してループ
            for (int k = sparse.row_pointers[j]; k < sparse.row_pointers[j + 1]; ++k) {
                int col_index = sparse.col_indices[k];
                double sparse_value = sparse.values[k];

                // 密行列のi行目のj列と疎行列の非ゼロ要素の掛け算
                temp[col_index] += dense.values[i][j] * sparse_value;
            }
        }

        // tempに蓄積された結果を疎行列に追加
        for (const auto& entry : temp) {
            result.addValue(i, entry.first, entry.second);
        }
    }

    result.build();  // 行ポインタを構築
    return result;
}

SparseMatrix operator*(const SparseMatrix& A, const SparseMatrix& B) {
    assert(A.n_col == B.n_row);  // 行列のサイズが適切であることを確認

    int rows = A.n_row;
    int cols = B.n_col;
    SparseMatrix result(rows, cols);

    // 一時的な結果を保持するハッシュマップ（列インデックスをキー、値が結果）
    std::unordered_map<int, double> temp;

    // 行列Aの各行についてループ
    for (int i = 0; i < A.n_row; ++i) {
        temp.clear();  // ハッシュマップをクリア

        // Aのi行目の非ゼロ要素に対してループ
        for (int k = A.row_pointers[i]; k < A.row_pointers[i + 1]; ++k) {
            int colA = A.col_indices[k];  // Aの非ゼロ要素の列インデックス
            double valueA = A.values[k];  // Aの非ゼロ要素の値

            // BのcolA行目の非ゼロ要素に対してループ
            for (int j = B.row_pointers[colA]; j < B.row_pointers[colA + 1]; ++j) {
                int colB = B.col_indices[j];  // Bの非ゼロ要素の列インデックス
                double valueB = B.values[j];  // Bの非ゼロ要素の値

                // 結果をハッシュマップに蓄積
                temp[colB] += valueA * valueB;
            }
        }

        // tempに蓄積された結果を疎行列に追加
        for (const auto& entry : temp) {
            if (entry.second != 0.0) {
                result.addValue(i, entry.first, entry.second);
            }
        }
    }

    result.build();  // 行ポインタを構築
    return result;
}

SparseMatrix operator * (double c, const SparseMatrix &x){
    SparseMatrix z = x;
    z *= c;
    return z;
}

Vector operator * (const Matrix &mat, const Vector &vec){ 
    assert(mat.n_col == vec.n_size);
    Vector y = Vector(mat.n_row, 0);
    for (size_t i = 0; i < mat.n_row; i++)
        for (size_t j = 0; j < mat.n_col; j++)
            y.values[i] += mat.values[i][j] * vec.values[j];
    return y;    
}

SparseVector operator * (const Matrix &mat, const SparseVector &vec){ 
    for(const auto [indice, value] : vec) assert(indice < mat.n_col);
    SparseVector y;
    for (size_t i = 0; i < mat.n_row; i++) { 
        double element = 0.0;
        for(const auto [indice, value] : vec)
            element += mat.values[i][indice] * value;
        if(element != 0.0)
            y.push(i, element);
    }
    return y;    
}

Matrix operator / (const Matrix &x, double c){
    Matrix z = x;
    z /= c;
    return z;
}

SparseMatrix operator / (const SparseMatrix &x, double c){ 
    SparseMatrix z = x;
    z /= c;
    return z;
}

std::ostream &operator << (std::ostream &stream, const Matrix &mat){
    
    for (size_t i = 0; i < mat.n_row; i++)
    {
        for (size_t j = 0; j < mat.n_col; j++)
          stream << mat.values[i][j] << ", ";
        stream << std::endl;
    }
    
    return stream;
}

std::ostream &operator << (std::ostream &stream, const SparseMatrix &mat){
    vec2d values(mat.n_row, std::vector<double>(mat.n_col, 0.0));
    for (size_t i = 0; i < mat.n_row; i++)
        for (size_t j = mat.row_pointers[i]; j < mat.row_pointers[i + 1]; ++j)
            values[i][mat.col_indices[j]] = mat.values[j];
    Matrix dense = Matrix(mat.n_row, mat.n_col, values);
    stream << dense;
    return stream;
}

// 非ゼロ要素を追加
void SparseMatrix::addValue(int row, int col, double value) {
    assert(row >= 0 && row < n_row);
    assert(col >= 0 && col < n_col);
    if (value != 0.0) {
        values.push_back(value);
        col_indices.push_back(col);
        row_pointers[row + 1]++;
    }
}

// 行ポインタを構築
void SparseMatrix::build() {
    for (size_t i = 1; i <= n_row; ++i) {
        row_pointers[i] += row_pointers[i - 1];  // 各行の累積
    }
}

// 要素を取得
double SparseMatrix::getValue(int row, int col) const {
    assert(row >= 0 && row < n_row);
    assert(col >= 0 && col < n_col);

    // 指定された行における非ゼロ要素を探索
    for (size_t i = row_pointers[row]; i < row_pointers[row + 1]; ++i) {
        if (col_indices[i] == col) {
            return values[i];
        }
    }
    return 0.0;  // 非ゼロ要素がない場合は0を返す
}